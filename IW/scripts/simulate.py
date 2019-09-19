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

#Source local libraries
sys.path.append('../src')

from iw_field import InternalWaveField
from iw_sim   import InternalWaveSimulation


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


#Frequency Distrubution (non radial)

freqs = np.array(p['freqs'])/cph
modes = np.array(p['modes'])
amps_real  = p['amps_real'] 
amps_imag  = p['amps_imag']
headings   = p['headings']

if len(amps_real) != len(modes)*len(freqs) != len(headings):
    print("Config file error: length amps != length modes*freqs")

amps = []
for i,a in enumerate(amps_real):
    zz = list(zip(a,amps_imag[i],headings[i]))
    amps.append(    { 'amps' : [ complex(z[0],z[1]) for z in zz],
                      'headings': [np.pi*z[2]/180 for z in zz]} )

#Make wave field
iwf = InternalWaveField(iwrange,iwdepth,
                        freqs=freqs,
                        modes=modes,
                        amplitudes=amps)


#Print parameters to users
print("INPUT PARAMETERS:")
print("\tFrequencies : " , freqs)
print("\tAmplitudes ",amps)
print("\tHeadings "  ,headings)
print("\tHorizontal Waves :", [iwf.hwavenumbers[m] for m in iwf.modes])

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
                
                
                
                

