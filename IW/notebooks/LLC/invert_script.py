#!/usr/bin/env python
# coding: utf-8



import os
import sys
import feather
import seawater as sw
from scipy import interpolate as interp
from scipy import signal as sig
sys.path.append("../scripts")
sys.path.append("../src/iw_model")
sys.path.append("../src/iw_inverse")
sys.path.append("../src/misc")

from iw_invert import InternalWaveInversion
from iw_plots import *
from extract import *

def process_data():
    """
    Read data and compute means
    """
    path = os.path.join("..","matlab",filename)
    table = read_table(path)
    
    #Define a slice of data
    DSL = slice(15,80) # depth slice
    sl = (slice(0,2000),DSL,slice(0,-1),slice(0,-1))
    
    #Get state variables
    N2 = extract_node(table,'N2')[sl]
    T  = extract_node(table,'T')[sl]
    S  = extract_node(table,'S')[sl]
    W  = extract_node(table,'W')[sl]*100
    Z  = extract_node(table,'z').flatten()[sl[1]]
    
    #Compute Density from state variables
    RP = sw.eos80.dens(S,T,2000)
    
    TIME   = extract_node(table,'time').flatten()[sl[0]]
    DEPTH  = extract_node(table,'z').flatten()[sl[1]]
    
    #Start time at zero
    t0 = matlab2datetime(TIME[1]).timestamp()
    TIME = [round(matlab2datetime(T).timestamp()-t0) for T in TIME]
    X = 1e3*np.arange(0,22,2)
    Y = 1e3*np.arange(0,22,2)
    
    #Form a coordinate matrix
    AX = np.array(list(itertools.product(*[TIME,DEPTH,X,Y])),dtype=[('t',float),('z',float),('x',float),('y',float)])
    AX = AX.reshape(T.shape)
    
    
    
    """
    Compute means
    """
    #Compute the mean stratification
    N2_mean = np.mean(N2,axis=(0,2,3))
    T_mean  = np.mean(T,axis=(0,2,3))
    S_mean  = np.mean(S,axis=(0,2,3))
    R_mean  = sw.eos80.dens(S_mean,T_mean,2000)
    #R_mean  = sw.eos80.pden(S_mean[DSL],T_mean[DSL],Z[DSL],pr=2000)
    
    
    #Create a function that is depth as a function of the mean density
    zofr = interp.InterpolatedUnivariateSpline(R_mean,Z)
    DZ   = AX['z'] - zofr(RP) 
    
    
    """
    M2 band pass filter
    """
    def bandpass_filter(center,half_width,fs,order=5):
        nyq = fs*0.5
        lo = (center-half_width/nyq)
        hi = (center+half_width/nyq)
        b, a = sig.butter(order, [lo, hi], btype='band')
        return b,a
    
    M2  = .0805 
    bw  = .02
    b,a = bandpass_filter(M2,bw,1)
    DZB = sig.lfilter(b,a,DZ,axis=0)
    RPB = sig.lfilter(b,a,RP,axis=0)
    WB = sig.lfilter(b,a,W,axis=0)
 
    return AX,DZ,DZB,Z,N2_mean

def compute_inversion(name,isl,AX,DZ,N2_mean):
    
    
    """
    Compute Inversion All Data
    """
    DSL = slice(15,80) # depth slice
    
    #Interpolate N2 to a regularized grid
    n2ofz = interp.InterpolatedUnivariateSpline(Z[DSL],N2_mean[DSL])
    Z_even = np.linspace(Z[10],Z[60],100)
    
    #Subsample data set create dataframe
    cube = pd.DataFrame({"dz"    : DZ[isl].flatten(),
                         "time"  : AX[isl]["t"].flatten(),
                         "z"     : AX[isl]["z"].flatten(),
                         "x"     : AX[isl]["x"].flatten(),
                         "y"     : AX[isl]["y"].flatten()
                        })
    
    
    #Compute depth variance then apply it to the dataset
    R = np.var(DZ[isl], axis=(0,2,3))
    Zr = AX[0,isl[1],0,0]["z"]
    snr = 0.9
    R = np.array([ (1-snr)*r if (1-snr)*r > 4 else 4 for r in R])
    RofZ = interp.interp1d(Zr,R)
    cube['RINV'] = 1/RofZ(cube['z'])
    
    #Set up parameter space
    FREQS = [.0805/3600]
    MODES = np.arange(0,5)
    ANGLES  = np.arange(0,360,5)
    
    #Run inversion
    subdf = cube[cube['time'] < 0.5*max(cube['time'])]
    iwi = InternalWaveInversion(subdf,FREQS,MODES,ANGLES,Z_even,n2ofz(Z_even),tc='dz')
    iwi.ps['QINV'] =(iwi.ps['modes'] + 3)**2# np.ones(len(iwi.ps))
    iwi.pinvert()
    
    print(name + " inversion is complete") 
    
    #Compute an estimate
    iwi_full = InternalWaveInversion(cube,FREQS,MODES,ANGLES,Z_even,n2ofz(Z_even),tc='dz')
    
    cube['dz_hat'] = iwi_full.H @ np.concatenate([iwi.ps['a'],iwi.ps['b']]).T
    cube['err']    = cube['dz_hat'] - cube['dz']
    
    print(name + " estimate  is complete") 
    
    
    #Write Out Dataframes
    iwi.ps =  iwi.ps.astype({'hwvn' : float})
    feather.write_dataframe(iwi.ps,"data/" + name+"_ps.fthr",)
    feather.write_dataframe(cube,  "data/" + name+"_ds.fthr")
    
    print(name + " program  is complete") 
    
"""
Main program
"""
if __name__== "__main__":
    AX,DZ,DZB,Z,N2_mean = process_data()
    
    """
    South axis test
    """
    isl = (slice(0,336,1),slice(10,50,2),slice(0,-1),slice(0,1))
    name = "south"
    compute_inversion(name,isl,AX,DZ,N2_mean)
    
    """
    West Y axis test
    """
    isl = (slice(0,336,1),slice(10,50,2),slice(0,1),slice(0,-1))
    name = "west"
    compute_inversion(name,isl,AX,DZ,N2_mean)
    
    """
    East X axis test
    """
    isl = (slice(0,336,1),slice(10,50,2),slice(0,-1),slice(-2,-1))
    name = "east"
    compute_inversion(name,isl,AX,DZ,N2_mean)
    
    """
    North Y axis test
    """
    isl = (slice(0,336,1),slice(10,50,2),slice(-2,-1),slice(0,-1))
    name = "north"
    compute_inversion(name,isl,AX,DZ,N2_mean)
    
    
    """
    Southwest corner
    """
    isl = (slice(0,336,1),slice(10,50,2),slice(0,4),slice(0,4))
    name = "swcorner"
    compute_inversion(name,isl,AX,DZ,N2_mean)
    
    
    """
    Northeast corner
    """
    isl = (slice(0,336,4),slice(10,50,2),slice(-5,-1),slice(-5,-1))
    name = "necorner"
    compute_inversion(name,isl,AX,DZ,N2_mean)

