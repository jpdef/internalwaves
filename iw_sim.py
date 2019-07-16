#IW_SIM
#Desc : Library for simulating internal wave fields
#Auth : J. DeFilippis
#Date : 7-16-2019

import feather
import os
from iw_field import InternalWaveField 
import numpy as np
import sys
from tqdm import tqdm


class InternalWaveSimulation:
    """
    Desc: 
    A class to generate time stepped interwave
    fields using numerical methods

    """

    def __init__(self,timeaxis,iwf,ftype=0,dpath="",fname=""):
        self.frames = []
        self.timeaxis = timeaxis
        self.iwf = iwf
        self.ftype = ftype
        self.dpath = dpath if dpath else os.getcwd()
        if not os.path.exists(self.dpath):
            os.mkdir(self.dpath)
        print(self.dpath)
        self.fname = fname if fname else "iwfsim"


    def run(self):
        print("Running simulation")
        for i,t in tqdm(enumerate(self.timeaxis),ascii=True,total=len(self.timeaxis),leave=True):
            step = self.make_step()
            self.iwf.update_field(step)
            self.frames.append(self.iwf.to_dataframe())
        sys.stdout.flush()
        print("\n")
        
        print("Writing to disk")
        self.make_files()
        sys.stdout.flush()
        print("\n")

    def make_step(self):
        waves = np.array( [np.exp(2*np.pi*1j*f) for f in self.iwf.freqs[0]])
        step = np.multiply(self.iwf.weights,waves)
        return step
    
    
    def make_files(self):
        if  self.ftype==0:
            self.make_featherfiles()
        
        elif self.ftype==1:
            self.make_csvfiles()
        
        elif self.ftype==2:
            self.make_animation()

    def make_featherfiles(self):
        for t,f in tqdm(enumerate(self.frames),ascii=True,total=len(self.frames),leave=True):
            fname = "%s-%d.fthr" % ( self.fname, t)
            path = os.path.join(self.dpath,fname)
            feather.write_dataframe(f,path) 


    def make_csvfiles(self):
        pass 

    
    def make_animation(self):
        pass


    def compute_file_size(self):
        pass


    def compute_run_time(self):
        pass
