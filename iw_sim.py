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
    A class to generate time stepped interwave
    fields using numerical methods

    """

    def __init__(self,timeaxis,iwf,ftype=0,dpath="",fname=""):
        self.frames = []
        self.fields = []
        self.scalarfields = []
        self.timeaxis = timeaxis
        self.iwf = iwf
        self.ftype = ftype
        self.delta_t = max(self.timeaxis)/( (len(self.timeaxis)-1) * (3600) )
        self.dpath = dpath if dpath else os.getcwd()
        self.zero_padding = int(np.floor( np.log10(len(self.timeaxis)) ) + 1)
        if not os.path.exists(self.dpath):
            os.mkdir(self.dpath)
        print("Datafile directory: ",self.dpath)
        self.fname = fname if fname else "iwfsim"


    def run(self):
        if len(self.timeaxis) > 1000:
            chunk_size = int( np.floor(len(self.timeaxis)/1000) )
            timechunks = self.make_chunks(chunk_size)
            for i,tc in self.progressbar(timechunks, "Long Simulation"):
                self.timeaxis = tc
                self.simulate()
                self.make_files(offset=i*chunk_size)
                self.fields = []
                self.frames = [] 
       
        else:
            self.simulate()
            self.make_files()
     
    def simulate(self):
        for i,t in self.progressbar(self.timeaxis,"Simulating"):
            step = self.make_step(t)
            self.iwf.update_field(step)
            self.frames.append(self.iwf.to_dataframe())
            self.fields.append(self.iwf.field)
            #self.scalarfields.append(self.iwf.scalarfield.field)


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
        waves = np.array( [np.exp(-2*np.pi*1j*f*t) for f in self.iwf.freqs[0]])
        step  = np.multiply(self.iwf.weights,waves)
        return step
    
    
    def make_files(self,offset=0):
        """
        Desc:
        Multiplexes throug hthe various output types
        """
        if  self.ftype==0:
            self.make_featherfiles(offset)
        
        elif self.ftype==1:
            self.make_csvfiles()
        
        elif self.ftype==2:
            self.make_animation()


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
        f = pd.DataFrame({'time'  : self.timeaxis,
                         'range' : self.iwf.range,
                         'depth' : self.iwf.depth,
                         'bfrq'  : self.iwf.bfrq})
        feather.write_dataframe(f,path) 


    def make_csvfiles(self):
        pass 

    
    def make_animation(self):
        fig, ax = plt.subplots()
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        dist  = self.iwf.range/1000
        depth = self.iwf.depth/1000
        self.set_animation_attributes(fig,ax) 
        
        #First Frame
        zeros  = np.zeros(shape=self.iwf.field.shape)
        p = ax.contourf(dist,depth,zeros,20,cmap=cmocean.cm.thermal)
        
        #Init and update function
        init = functools.partial(self.init_animation,fig,ax,p,cbar_ax)
        update = functools.partial(self.update_animation,fig,ax,cbar_ax)
       
        #Compile Animation
        print("Creating Animation it may take awhile")
        ani = FuncAnimation(fig, update, frames=np.arange(0,len(self.timeaxis),1),init_func=init, blit=False)
        path = os.path.join(self.dpath,'animation.mp4')
        ani.save(path)


    def set_animation_attributes(self,fig,ax):
        ax.invert_yaxis()
        ax.set_xlabel('Range Km')
        ax.set_ylabel('Depth Km')
        
    def init_animation(self,fig,ax,p,cbar_ax):
        fig.subplots_adjust(right=0.8)
        cbar = fig.colorbar(p, cax=cbar_ax)
        cbar.set_label('Displacement (m)')
        return p

    def update_animation(self,fig,ax,cbar_ax,frame):
        ax.set_title('%.2f hours' % np.multiply(frame,self.delta_t))
        field = self.fields[frame]
        #field = self.scalarfields[frame].real
        dist  = self.iwf.range
        depth = self.iwf.depth
        p = ax.contourf(dist,depth,field.real,20,cmap=cmocean.cm.thermal)
        cbar_ax.cla()
        fig.colorbar(p, cax=cbar_ax)
        return p

    def compute_file_size(self):
        pass


    def compute_run_time(self):
        pass
