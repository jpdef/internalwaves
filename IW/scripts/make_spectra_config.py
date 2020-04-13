#!/usr/bin/python3
#SIMULATE_SPECTRA
#Desc : Script that executes internal wave simulation
#       based on spectral parameters
#Auth : J. DeFilippis
#Date : 2-29-2020

import sys
import numpy as np
import pandas as pd
import feather
import json
import os 

#Source local libraries
sys.path.append('../src/iw_model')
sys.path.append('../src/misc')

from iw_modes import InternalWaveModes


def psd(f,sigma_sq):
    return sigma_sq - 2*f/max(f)
    
    
def prompt_config():
    #Get config file from user
    if len(sys.argv) > 1:
        config_fname = sys.argv[1]
    else:
        print("Usage : need configuration filename ")
        sys.exit(1)
    
    #Read Sim Params
    with open(config_fname) as param_file:
        p = json.load(param_file)
    
    return p

p = prompt_config()


"""
Frequency range from fundmental (1/T) to Nyquist
   -If fundamental is below interia frequency assign default to 23 hour
   -Longer time simulation -> better spectra resolution
   -Smaller time step more frequencies in spectra
"""
try:
    freqs = np.array( p['freqs'] )/3600
except KeyError:
    T = (3600*23)  if (p['time_stop'] > 3600*24)  else p['time_stop']
    ff = 1/T
    fn = 1/(2*p['time_step']) if( 1/(2*p['time_step']) < (4/3600) ) else (4/3600)
    freqs = np.arange(ff,fn,ff)


"""
Wavenumber range is determined from frequency range using the
frequency range
    -Compute wave directions from dk resolution from nyquist frequency
"""
depth = np.linspace(0,p['depth_end'],p['depth_res'])
N2    = feather.read_dataframe(p['envfile'])['strat']
K = []
for i in range(len(freqs)):
    iwm = InternalWaveModes(depth,N2,freq=freqs[i])
    K.append([iwm.get_hwavenumber(m) for m in p['modes']] )

K = np.array(K).flatten()
print(K)
p['K'] = list(K.real)


#Fundamental wavenumber
dk = K[0].real

headings = []
for kmag in K:
    dtheta = 360 * (dk/kmag.real)
    headings.append( np.arange(0,360,dtheta) )
headings = [list(h) for h in headings]


"""
Update amplitudes for headings
"""

#Compute amplitude and duplicate for each mode1
amps  = psd(freqs,1e-2)
amps  = np.repeat(amps,len(p['modes']))

#Random the phase
phi = np.random.uniform(low=0,high=2*np.pi,size=len(amps))
amps_real = np.multiply(amps,np.sin(phi))
amps_imag = np.multiply(amps,np.cos(phi))

#Normalize for the amount of directions
amps_real =  [ np.repeat(a/sqrt(len(headings[i])) ,len(headings[i])) for i,a in enumerate(amps_real) ]
amps_imag =  [ np.repeat(a/sqrt(len(headings[i])) ,len(headings[i])) for i,a in enumerate(amps_imag) ]


"""
Print out Diagnostics and suggested sampling for inversion
"""

p['freqs']     = list(3600*freqs)
p['amps_real'] = [list(a) for a in amps_real]
p['amps_imag'] = [list(a) for a in amps_imag]
p['headings']  = headings


with open('out.json','w') as out_file:
    json.dump(p,out_file)
    
