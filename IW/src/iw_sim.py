#IW_SIM
#Desc : Library for simulating internal wave fields
#Auth : J. DeFilippis
#Date : 7-16-2019

import feather
import os
import pandas as pd
from iw_field import InternalWaveField 
import numpy as np
import sys
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmocean
import functools

class InternalWaveSimulation:
    """
    Desc: 
    A class to generate time stepped internalwave
    fields using numerical methods

    """

    def __init__(self,timeaxis,iwf,ftype=0,dpath="",fname="",chunklim=100):
        self.frames = []
        self.fields = []
        self.timeaxis = timeaxis
        self.iwf = iwf
        self.ftype = ftype
        self.chunklim = chunklim
        self.delta_t = max(self.timeaxis)/( (len(self.timeaxis)-1) * (3600) )
        self.dpath = dpath if dpath else os.getcwd()
        self.zero_padding = int(np.floor( np.log10(len(self.timeaxis)) ) + 1)
        if not os.path.exists(self.dpath):
            os.mkdir(self.dpath)
        print("Datafile directory: ",self.dpath)
        self.fname = fname if fname else "iwfsim"


    def run(self,coords=[]):
        if len(self.timeaxis) > self.chunklim and self.ftype != 2:
            chunk_size = int( np.floor(len(self.timeaxis)/self.chunklim) )
            timechunks = self.make_chunks(chunk_size)
            for i,tc in self.progressbar(timechunks, "Long Simulation"):
                self.timeaxis = tc
                self.simulate(coords=coords)
                self.make_files(offset=i*chunk_size)
                self.frames = []
                self.fields=[] 
       
        else:
            self.simulate(coords=coords)
            self.make_files()
     
    def simulate(self,coords=[]):
        for i,t in self.progressbar(self.timeaxis,"Simulating"):
            step = self.make_step(t)
            self.iwf.update_field(step)
            self.frames.append(self.iwf.to_dataframe(coords=coords,time=t))
            self.fields.append(self.iwf.field)


    def progressbar(self,dataset,desc):
        """
        Desc:
        Helper function that wraps the tqdm library to make 
        function call shorter
        """
        iterator = enumerate(dataset)
        return tqdm(iterator,ascii=True,total=len(dataset),leave=True,desc=desc)

 
    def make_chunks(self,chunk_size):
        timechunks = []
        back_itr = 0
        forward_itr=chunk_size
        while forward_itr < len(self.timeaxis):
             timechunks.append( self.timeaxis[back_itr : forward_itr] )
             back_itr=forward_itr
             forward_itr += chunk_size
        #Remainder
        timechunks.append( self.timeaxis[back_itr:] ) 

        return timechunks


    def make_step(self,t):
        """
        Desc:
        Modulates each frame with a e^(2pi i * f *t) where t is the
        time input and the frequencies are known.
        """
        waves = [] 
        for f in self.iwf.freqs:
            waves.append(np.exp(-2*np.pi*1j*f*t))
        waves = np.array(waves)
        return waves
    
    
    def make_files(self,offset=0):
        """
        Desc:
        Multiplexes throug hthe various output types
        """
        if  self.ftype==0:
            self.make_featherfiles(offset)
        
        elif self.ftype==1:
            self.make_csvfiles()
        

    def make_featherfiles(self,offset=0):
        for t,f in self.progressbar(self.frames,"Writing to Disk"):
            fmt = '{:0>' + str(self.zero_padding) + '}'
            fname = "%s-%s.fthr" % ( self.fname, fmt.format(t+offset) )
            path = os.path.join(self.dpath,fname)
            feather.write_dataframe(f,path) 


    def make_metadata_file(self):
        #Write a file contain various parameters from simulation
        fname = "meta.fthr"
        path = os.path.join(self.dpath,fname)
        f = pd.DataFrame({'time_len'  : float(len(self.timeaxis)),
                          'time_max'  : float(max(self.timeaxis)),
                          'range_len' : float(len(self.iwf.range)),
                          'range_max' : float(max(self.iwf.range)),
                          'depth_len' : float(len(self.iwf.depth)),
                          'depth_max' : float(max(self.iwf.depth))},
                         index=[0])
                         #'bfrq'  : self.iwf.bfrq},
        feather.write_dataframe(f,path) 


    def make_csvfiles(self):
        pass 

    

    def compute_file_size(self):
        pass


    def compute_run_time(self):
        pass
