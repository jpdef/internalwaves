#!/usr/bin/env python


import scipy.io as scio
import pandas as pd 
import xarray as xr
import datetime as dt
import utils
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import rc
import scipy.signal as signal
import plotly
import numpy as np

def convert_to_ts(seconds):
    dt_list = []
    for s in np.nditer(seconds):
        dt_list.append(datetime.fromtimestamp(s))
    return dt_list

def load_dataframe(path):
    keys = ['time','depth','temp']
    mat  = scio.loadmat(path)
    return pd.DataFrame(dict((k, mat[k].flatten()) for k in keys))

def plot_psd_(ax,D,M,FS,label,si=slice(1,None)):
    freqs,psd   = utils.psd(D,M=M,FS=FS)
    high,low    = utils.psd_err(M) #Specific for welch with hanning
    ax.loglog(freqs[si],psd[si],label=label)
    ax.fill_between(freqs[si],low*psd[si],high*psd[si],alpha=0.2)
    return freqs,psd


def plot_coh(ax,D1,D2,M,FS,label,si=slice(1,None)):
    nperseg = int(len(D1)/M)
    nd = 2*M*0.95
    freqs,coh = signal.coherence(D1,D2,nperseg=nperseg,fs=FS)
    beta = 1 - (0.005)**(1/(nd-1))
    ax.semilogx(freqs[si],coh[si])
    ax.plot(freqs[si],beta*np.ones(len(freqs[si])),'r--',label='Sig. Threshold 99')
    ax.set_xlabel('Frequency (CPD)')
    ax.set_ylabel(r'$\gamma^2$')
    ax.grid()
    ax.legend()
    return freqs[si],coh[si],beta


def plot_csd(ax,D1,D2,M,FS,label,si=slice(1,None)):
    nperseg = int(len(D1)/M)
    freqs, Pxy = signal.csd(D1,D2, nperseg=nperseg,fs=FS)
    ax.plot(freqs[si],Pxy.real[si],marker='s',label='real')
    ax.plot(freqs[si],Pxy.imag[si],marker='s',label='imag')
    ax.grid()
    ax.legend()
    return freqs[si],Pxy[si]

def plot_phase(ax,D1,D2,M,FS,label,si=slice(1,None)):
    nperseg = int(len(D1)/M)
    freqs,coh  = signal.coherence(D1,D2,nperseg=nperseg,fs=FS)
    freqs, Pxy = signal.csd(D1,D2, nperseg=nperseg,fs=FS)
    phase = np.angle(Pxy[si],deg=False)
    err = np.sqrt(1-coh[si]/(2*M))/(np.sqrt(coh[si]))
    ax.errorbar(freqs[si],phase,yerr=err,marker='s',ls='none',label=label)
    ax.set_ylabel('Radians')
    ax.grid()
    ax.legend()
    return freqs[si],phase,err


#Thermistor chain at 34.2545, -120.0974
TS1 = load_dataframe('../matlab/TS_1.mat')

#Thermistor chain at 34.2389, -120.1026
TS2 = load_dataframe('../matlab/TS_2.mat')

#Depth Average
TS1_DA = TS1.groupby('time').mean()
TS2_DA = TS2.groupby('time').mean()

#Spectral Parameters
FS = 60*24 #Sampled per minute meaning 24*60 samples per day
M  = 4      #Number of segments

#Plot the Coh,CSD,Phase
f,ax = plt.subplots(3,1)
freqs,coh,beta  = plot_coh(  ax[0],TS1_DA.temp,TS2_DA.temp,M,FS,'Sq. Coh')
freqs,csd       = plot_csd(  ax[2],TS1_DA.temp,TS2_DA.temp,M,FS,'CSD',slice(1,10))
freqs,phase,err = plot_phase(ax[1],TS1_DA.temp,TS2_DA.temp,M,FS,'Phase',slice(1,10))

#Title
f.suptitle('Square Coherence')
f.set_size_inches(10,7)
f.set_dpi(200)


sig_ind = np.where(coh>beta)[0][1:10]
f,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)
freqs,coh,beta  = plot_coh(  ax1,TS1_DA.temp,TS2_DA.temp,M,FS,'Sq. Coh',sig_ind)
freqs,csd       = plot_csd(  ax2,TS1_DA.temp,TS2_DA.temp,M,FS,'CSD',sig_ind)
freqs,phase,err = plot_phase(ax3,TS1_DA.temp,TS2_DA.temp,M,FS,'Phase',sig_ind)

wl = 4*np.pi/abs(phase)
wl_err = abs(wl)*(err/abs(phase))
ax4.plot(freqs, wl,marker='s',ls='none')


ax4.set_xlabel('Frequency (CPD)')
ax4.set_ylabel(r'Wavelength $\lambda$')
ax4.grid(which='both')

f.suptitle('Coherent Frequency and Phase')
f.set_size_inches(10,7)
f.set_dpi(200)
plt.show()
#
#
print(phase*180/np.pi)
print(freqs/24)


