#!/usr/bin/python3
#SIMULATE
#Desc : Script that executes internal wave simulation 
#Auth : J. DeFilippis
#Date : 7-16-2019

import sys
import numpy as np
import pandas as pd
import feather
import json
import os 
import functools
import scipy.interpolate as interp

#Source local libraries
sys.path.append('../src/iw_model')
sys.path.append('../src/misc')

from iw_field import InternalWaveField
from iw_sim   import InternalWaveSimulation
from map_scalars import map_scalars
from map_scalars import map_sound_speed
from map_scalars import map_sound_speed_anom

cph = 3600

#Get config file from user
if len(sys.argv) > 1:
    config_fname = sys.argv[1]
else:
    print("Usage : need configuration filename ")
    sys.exit(1)


#Read Sim Params
with open(config_fname) as param_file:
    p = json.load(param_file)

#Spacial Params
iwrange = np.linspace(0,p['range_end'],p['range_res'])
iwdepth = np.linspace(0,p['depth_end'],p['depth_res'])

#Stratification Profile
strat=np.array([])
try :
    df = feather.read_dataframe(p['n2_file'])
    strat = df['strat']
except KeyError:
    pass
     
#Frequency Distrubution (non radial)

freqs = np.array(p['freqs'])/cph
modes = np.array(p['modes'])
amps_real  = p['amps_real'] 
amps_imag  = p['amps_imag']
headings   = p['headings']

if len(amps_real) != len(modes)*len(freqs) != len(headings):
    print(len(amps_real), len(modes),len(freqs), len(headings))
    print("Config file error: length amps != length modes*freqs")

amps = []
for i,a in enumerate(amps_real):
    zz = list(zip(a,amps_imag[i],headings[i]))
    amps.append(    { 'amps' : [ complex(z[0],z[1]) for z in zz],
                      'headings': [2*np.pi*z[2]/360 for z in zz]} )


#Make wave field
iwf = InternalWaveField(iwrange,iwdepth,
                        freqs=freqs,
                        modes=modes,
                        amplitudes=amps,
                        bfrq=strat)


#Print parameters to users
print("INPUT PARAMETERS:")
print("\tFrequencies : " , freqs)
#print("\tHeading & Amplitude : ",amps)
print("\tHorizontal Wavenumber :", [iwf.iwmodes[0].get_hwavenumber(m) for m in iwf.modes])


#Make sampling coordinates
coords = []
for xi in p['x_samples']:
    for yi in p['y_samples']:
        for zi in p['z_samples']:
            coords.append( (xi,yi,zi) )


#Run simulation
time = np.arange(0,p['time_stop']+p['time_step'],p['time_step'])
iws = InternalWaveSimulation(time,iwf=iwf,dpath=p['path'],fname='run',ftype=p['ftype'])
iws.make_metadata_file()
iws.run(coords=coords) 
                
                
                
#Mapping Sound Profile onto IW Field
if p['ssp_file']:
    strat=np.array([])
    try :
        df = feather.read_dataframe(p['ssp_file'])
        ssp = df['ssp']
        ssz = df['depth']
    except KeyError:
        print("Incorrect ssp format")
        exit(0)
   
    #Compute interpolation and derivative
    # TODO needs to be potential sound speed gradient
    fssp  = interp.UnivariateSpline(ssz, ssp)
    fsspz = fssp.derivative()

    #include into partial function
    mss  = functools.partial(map_sound_speed,ssp=fssp,ssp_grad=fsspz)
    mssa = functools.partial(map_sound_speed_anom,ssp_grad=fsspz)
    map_scalars(p['path'],'c',mss)
    map_scalars(p['path'],'ca',mssa)


else:
    map_scalars(p['path'],'c',map_sound_speed)
    map_scalars(p['path'],'ca',map_sound_speed_anom)
