#Desc : Library for evaluating parameters of observed internal wave field
#Auth : J. DeFilippis
#Date : 10-28-2019

import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm

import sys
sys.path.append('../src')

import pseudo_inverse as pinv

from iw_modes import InternalWaveModes


class InternalWaveInversion:
    """
    Desc: Class to invert for parameters frequency,mode, and wavelength
    from oceanographic data using a linear internal wave model
    """
    
    def __init__(self,ds,freqs,modes,thetas,depth,
                 N2=np.array([]),hwvn=np.array([]),tc='dz',mean=True):
        """
        Desc: Sets inputs to field values. Creates a dataframe of parameters
        to invert for. Creates a model matrix from data coordinates x,y,z,t.
        """
        self.set_attributes(ds,freqs,modes,thetas,depth,N2,hwvn,tc)
        self.ps = self.generate_param_space()
        self.H  = self.make_model_matrix()
        
        #Fit a mean if the frequency is zero
        if mean:
            new_row = pd.DataFrame({"freqs" : [0], "modes" : [0], "hwvn" : [0], "theta" : [0]})
            self.ps = pd.concat([new_row, self.ps[:]], sort=True).reset_index(drop = True) 
            
            mid = int( np.floor((self.H.shape[1])/2) ) + 1
            self.H = np.insert(self.H,0,np.ones(self.H.shape[0]).T,axis=1)
            self.H = np.insert(self.H,mid,np.ones(self.H.shape[0]).T,axis=1)
    
    
    def set_attributes(self,ds,freqs,modes,thetas,depth,N2,hwvn,tc):
        self.freqs  = freqs
        self.modes  = modes
        self.N2     = N2
        self.depth  = depth
        self.ds     = ds
        self.tc     = tc
        self.thetas = thetas
        self.hwvn   = hwvn
    
    
    def generate_param_space(self):
        #Create combination of frequency,modes, and wave direction
        if len(self.hwvn) > 0 :
            parr = np.array(list(itertools.product(*[self.freqs,self.hwvn,self.modes,self.thetas])),
                         dtype=[('f',float),('k',float),('m',int),('theta',float)])
        else:
            parr = np.array(list(itertools.product(*[self.freqs,self.modes,self.thetas])),
                         dtype=[('f',float),('m',int),('theta',float)])
        #Construct data frame
        param_df = pd.DataFrame({'freqs' : parr['f'], 
                                 'modes' : parr['m'], 
                                 'theta' : parr['theta']
                            })
        
        #Generate vertical modal structures
        self.iwmodes = [InternalWaveModes(self.depth,self.N2,freq=f) for f in self.freqs]
       
       
        #Fill in wave numbers
        param_df['hwvn'] = parr['k'] if len(self.hwvn) > 0 else self.map_wavenumber(parr) 
        
        return param_df 
    
    
    def map_wavenumber(self,parr):    
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx
    
        #Map wavenumbers to mode,frequency
        hwvn = []
        for p in parr:
            fi = find_nearest(self.freqs,p['f'])
            hwvn.append(self.iwmodes[fi].get_hwavenumber(p['m'])) 

        return hwvn
    
            
    def make_model_matrix(self):
        vmat = vertical_matrix(self.ds,self.ps,self.iwmodes)
        hmat = horizontal_matrix(self.ds,self.ps)
        
        return np.multiply(vmat,hmat)
    

    def pinvert(self):
        if 'RINV' in self.ds:
            RINV = np.array(self.ds['RINV'])  
        else:
            RINV = None 
        
        if 'QINV' in self.ps:
            QINV = np.concatenate([self.ps['QINV'],self.ps['QINV']])
        else:
            QINV = None
        
        self.hdag = pinv.tappered_least_square(self.H,RINV=RINV,QINV=QINV)
            
        amps = self.hdag @ self.ds[self.tc] 
        div  = int(np.floor(len(amps)/2))
        self.ps['a'] = amps[0:div]
        self.ps['b'] = amps[div:]

        
    def pinvert_RLS(self,itr_ax_label='time',li=1):
        itr_ax = np.unique(self.ds[itr_ax_label])
        ab = np.random.rand(2*len(self.ps))
        Q = np.diag(np.concatenate([self.ps['Q'],self.ps['Q']]))

        for it in tqdm(itr_ax):
            #Grab time slice
            sub_ds = self.ds.loc[self.ds[itr_ax_label]==it]
            
            #Update H matrix & Predict
            H = model_matrix(sub_ds,self.ps,self.iwmodes)
            est = H @ ab
            innovation = sub_ds['dz'] - est
            
            #Correct & Update
            R = np.diag(sub_ds['R'].values)
            INV = np.linalg.inv(li * H @ Q @ H.T + R)
            K = Q @ H.T @ INV
            ab = ab + K @ innovation
            Q = li*(Q - K @ H @ Q)
        
        #Update final parameters
        l = len(self.ps)
        self.ps['a'] = ab[:l]
        self.ps['b'] = ab[l:]
    
    
    def update_params(self, ps):
        self.ps = pd.concat([self.ps,ps],sort=True)
        Hnew = model_matrix(self.ds,ps,self.iwmodes)
        self.H = np.concatenate([self.H,Hnew],axis=1)
    
        
    def update_data(self,ds):
        self.ds = pd.concat([self.ds,ds],sort=True)
        Hnew = model_matrix(ds,self.ps,self.iwmodes)
        self.H = np.concatenate([self.H,Hnew],axis=0)


    def estimate(self,ds=None):
        if ds is None:
            return self.H @ np.concatenate([self.ps['a'],self.ps['b']]).T
        else:
            Hnew = model_matrix(ds,self.ps,self.iwmodes)
            return Hnew @ np.concatenate([self.ps['a'],self.ps['b']]).T
        
"""
Plane Wave + Model Matrices
"""    
def model_matrix(ds,ps,iwmodes):
    vmat = vertical_matrix(ds,ps,iwmodes)
    hmat = horizontal_matrix(ds,ps)
    
    return np.multiply(vmat,hmat)
    
    
def horizontal_matrix(ds,ps):
    """
    Desc: Create horizontal portion of matrix with sin and cos basis functions. 
    ie :// a sin(k . x  - ft) + b cos( k . x - ft)
    """
    #Make arguement
    rad =  ps['theta'] * (np.pi/180)
    kx = 2*np.pi*np.outer(ds['x'],ps['hwvn'].values.real*np.cos(rad)) 
    ky = 2*np.pi*np.outer(ds['y'],ps['hwvn'].values.real*np.sin(rad)) 
    ft = 2*np.pi*np.outer(ds['time'], ps['freqs'].values.real)
    
    #Compute for sin and cos parts
    psi  = np.hstack([np.sin(kx + ky - ft),np.cos(kx + ky - ft)])
    psi = psi.astype('float')
    return psi


def vertical_matrix(ds,ps,iwmodes):
    """
    Desc: Create  portion of matrix with vertical basis functions. 
    ie :// phi_nm (z)
    """
    cols = []
    
    #Get index for frequency values
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx
    
    #Iterate over parameter rows to get correct modal function
    #put them in a list
    freqs = [iwmode.freq for iwmode in iwmodes]
    for index,row in ps.iterrows():
        f = row['freqs'] 
        m = int(row['modes'].real)
        i = find_nearest(freqs,f)
        col = iwmodes[i].evaluate_mode(m,ds['z'])
        cols.append(col)
        
    """"
    Make a matrix rows of phi(z) then duplicate for cosine part
    and transpose so parameters are in column space and data
     is in row space 
    """
    phi = np.vstack(cols)
    phi = np.vstack([phi,phi])
    phi = phi.T
    phi = phi.astype(float)
    return phi