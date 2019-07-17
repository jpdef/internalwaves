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
        self.timeaxis = timeaxis
        self.iwf = iwf
        self.ftype = ftype
        self.delta_t = max(self.timeaxis)/(len(self.timeaxis)-1)
        self.dpath = dpath if dpath else os.getcwd()
        if not os.path.exists(self.dpath):
            os.mkdir(self.dpath)
        print(self.dpath)
        self.fname = fname if fname else "iwfsim"


    def run(self):
        print("Running simulation")
        for i,t in tqdm(enumerate(self.timeaxis),ascii=True,total=len(self.timeaxis),leave=True):
            step = self.make_step(t)
            self.iwf.update_field(step)
            self.frames.append(self.iwf.to_dataframe())
            self.fields.append(self.iwf.field)
        sys.stdout.flush()
        print("\n")
        
        print("Writing to disk")
        self.make_files()
        sys.stdout.flush()
        print("\n")

    def make_step(self,t):
        waves = np.array( [np.exp(-2*np.pi*1j*f*t) for f in self.iwf.freqs[0]])
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
            zero_padding = int(np.floor( np.log10(len(self.timeaxis)) ) + 1)
            fmt = '{:0>' + str(zero_padding) + '}'
            fname = "%s-%s.fthr" % ( self.fname, fmt.format(t))
            path = os.path.join(self.dpath,fname)
            feather.write_dataframe(f,path) 


    def make_csvfiles(self):
        pass 

    
    def make_animation(self):
        fig, ax = plt.subplots()
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        dist  = 10*self.iwf.range
        depth = 10*self.iwf.depth
        self.set_animation_attributes(fig,ax) 
        
        #First Frame
        zeros  = np.zeros(shape=self.iwf.field.shape)
        p = ax.contourf(dist,depth,zeros,20,cmap=cmocean.cm.thermal)
        
        #Init and update function
        init = functools.partial(self.init_animation,fig,ax,p,cbar_ax)
        update = functools.partial(self.update_animation,fig,ax,cbar_ax)
       
        #Compile Animation
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
        dist  = 10*self.iwf.range
        depth = 10*self.iwf.depth
        p = ax.contourf(dist,depth,field.real,20,cmap=cmocean.cm.thermal)
        cbar_ax.cla()
        fig.colorbar(p, cax=cbar_ax)
        return p

    def compute_file_size(self):
        pass


    def compute_run_time(self):
        pass
