#Desc : Library for computing internal wave fields
#Auth : J. DeFilippis
#Date : 6-26-2019

import matplotlib.pyplot as plt
from iw_modes  import InternalWaveModes
from iw_param  import Output
import scipy
import feather
import pandas as pd
from scipy import interpolate
import numpy as np

class InternalWaveField:
    """
    Desc: 
    A class to generate internal wave fields
    """

    def __init__(self,iwrange,iwdepth,
                 freqs=np.array([]), 
                 modes=np.array([0]),
                 amplitudes=[],
                 bfrq=np.array([]), 
                 f=1.1583e-5,
                 offset=[0,0]):
        
        print("Intializing wavefield")
        
        #Set all fields in object
        self.set_attributes(bfrq,iwrange,iwdepth,modes,f,offset)        
        
        #Compute phase speeds and vertical structure functions
        self.init_dispersion(freqs,amplitudes)

        #Compute 3D field from phase speeds and structure functions
        self.field = self.construct_field()
        
    
    def construct_field(self,step=np.array([])):
        """
        Desc:
        Constructs a 3D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        field = self.empty_field()
                
        #Initialize Field Components 
        if not self.field_components:
            for n in range(self.nfreqs):
                self.field_components.append(
                     self.construct_field_component(n))
      
        #Update Field Component Values
        for n,fc in enumerate(self.field_components):
           if step.size:
               field['w'] += fc['w']*step[n] 
               field['u'] += fc['u']*step[n] 
               field['v'] += fc['v']*step[n] 
               field['p'] += fc['p']*step[n] 
               field['d'] += fc['d']*step[n] 
           else:
               field['w'] += fc['w']
               field['u'] += fc['u']
               field['v'] += fc['v']
               field['p'] += fc['p']
               field['d'] += fc['d']
        
        return field   
    
    
    def construct_field_component(self,n):
        """
        Desc:
        Constructs a 3D wave field componenent for a specific frequency/wavenumber
        depending on the choosen dependent variable. These components are then summed
        to give the total wavefield.
        """
        field = self.empty_field()
       
        for nm,m in enumerate(self.modes):
            psi,psi_u,psi_v  = self.horizontal_comp(m,n)
            d    = self.iwmodes[n].d_modes[m]
            p    =  1j*self.iwmodes[n].p_modes[m]
            u    =  1j*self.iwmodes[n].u_modes[m]  
            w    = -1j*self.iwmodes[n].d_modes[m] * self.freqs[n]
            for zn in range(len(self.depth)):
                for xn in range(psi.shape[0]):
                    for yn in range(psi.shape[1]):
                        field[zn,xn,yn]['d']  += psi[xn,yn]*d[zn]
                        field[zn,xn,yn]['u']  += psi_u[xn,yn]*u[zn]
                        field[zn,xn,yn]['v']  += psi_v[xn,yn]*u[zn]
                        field[zn,xn,yn]['p']  += psi[xn,yn]*p[zn]
                        field[zn,xn,yn]['w']  += psi[xn,yn]*w[zn]  
      
        return field
   
    
    def horizontal_comp(self,nm,nf):
        """
        Desc:
        Constructs a plane wave in the horizontal with a specific 
        wavenumber and heading for mode index nm and frequency index nf 
        """
        xx,yy  = np.meshgrid(self.range - self.offset[0],self.range - self.offset[1])
        kmag   = self.iwmodes[nf].get_hwavenumber(nm)
        sqsum  = np.sqrt(self.freqs[nf]**2 + self.f**2)
        r1     = self.freqs[nf]/sqsum
        r2     = self.f/sqsum
        amps   = self.get_amplitude(nm,nf)
        psi    = np.zeros ( shape=xx.shape, dtype='complex128' )
        psi_u  = np.zeros ( shape=xx.shape, dtype='complex128' )
        psi_v  = np.zeros ( shape=xx.shape, dtype='complex128' )
        for ah in list(zip(amps['amps'],amps['headings'])):
            pw     = self.plane_wave(kmag,ah[0],ah[1],xx,yy)
            psi    += pw
            psi_u  += pw*(r1*np.cos(ah[1]) + 1j*r2*np.sin(ah[1])) 
            psi_v  += pw*(r2*np.sin(ah[1]) - 1j*r2*np.cos(ah[1]))
        
        return(psi,psi_u,psi_v)

  
    def plane_wave(self,kmag,amp,heading,xx,yy):
        kx      = kmag*np.cos(heading)
        ky      = kmag*np.sin(heading)
        psi     = amp*np.exp(2*np.pi*1j*(kx*xx +ky*yy)) 
        return(psi)
   
     
    def get_amplitude(self,m,n):
        """
        Desc:
        Returns amplitude for wave with mode m and frequency n
        """
        return self.amplitudes[m*self.nfreqs + n]


    def empty_field(self):
        """
        Desc:
        """ 
        field = np.zeros(shape=(len(self.range),len(self.range),len(self.depth))
                        ,dtype=[('d','complex'),('w','complex') ,('p','complex'),
                                ('u','complex'),('v','complex') ])
        return field


    def default_amplitudes(self,n):
        return([{'ar': 1 ,'ai' : 0 ,'theta' : np.pi/4} for n in range(n)])


    def init_dispersion(self,freqs,amplitudes):
        """
            z_n depth grid point for vertical functions
            m   mode number
            f_n associated frequency
        """
        D = len(self.depth)
        self.nwaves = len(freqs)*len(self.modes)
        self.nfreqs = len(freqs)
        
        
        self.freqs         = freqs
        self.iwmodes       = self.construct_modes(freqs)
        self.amplitudes    = amplitudes if amplitudes else self.default_amplitudes(self.nwaves)


    def construct_modes(self,freqs):
        """
        Desc:
        Constructs a set of modes from the set of input frequencies 
        via the dispersion relations. Mode objects contain wavenumbers
        """
        iwmodes = []
        print(freqs)
        for i in range(len(freqs)):
            iwmodes.append(InternalWaveModes(self.depth,self.bfrq,freq=freqs[i]))
        return iwmodes


    def update_field(self,step):
        """
        Desc:
        Efficient way to timestep or change distrubution after having calculated all
        the horizontal and vertical components
        """
        self.field = self.construct_field(step=step)

 
    def set_attributes(self,bfrq,iwrange,iwdepth,modes,f,offset):
        """
        Desc:
        Helper function for constructor to set all the various fields
        """
        self.range   = iwrange
        self.depth   = iwdepth
        self.zdiff   = np.average(np.diff(iwdepth)) 
        self.bfrq    = bfrq 
        self.modes   = modes
        self.f       = f
        self.offset  = offset
        self.field_components = [] 
    
    def to_dataframe(self,coords=[],time=0):
        """
        Desc:
        Converts the 3D array internal wave field to a panda's
        dataframe for use and convience
        """
    
        return self.select_data(coords,time) if coords else self.flatten_data(time) 
   
    
    def select_data(self,coords,time):
         x = [self.range[c[0]] for c in coords] 
         y = [self.range[c[1]] for c in coords] 
         z = [self.depth[c[2]] for c in coords]  
         d = [self.field[c[2],c[1],c[0]]['d'].real for c in coords]
         p = [self.field[c[2],c[1],c[0]]['p'].real for c in coords]
         u = [self.field[c[2],c[1],c[0]]['u'].real for c in coords]
         v = [self.field[c[2],c[1],c[0]]['v'].real for c in coords]
         w = [self.field[c[2],c[1],c[0]]['w'].real for c in coords]
         t = np.repeat(time,len(x))
        
         return pd.DataFrame({"x" : x , "y" :  y , "z" : z ,
                              "t" : t , "d" : d, "p" : p,
                              "u" : u , "v" : v, "w" : w})
    

    def flatten_data(self,time):
        x    = np.ndarray(shape=(1,),dtype=float)
        y    = np.ndarray(shape=(1,),dtype=float)
        z    = np.ndarray(shape=(1,),dtype=float)
        d = self.field['d'].real.flatten()
        p = self.field['p'].real.flatten()
        u = self.field['u'].real.flatten()
        v = self.field['v'].real.flatten()
        w = self.field['w'].real.flatten()
       
        L = len(self.range)
        H = len(self.depth)

        x = np.tile(self.range,L*H) 
        t = np.repeat(time,len(x))
           
        for i,yi in enumerate(self.range):
            y = yi*np.ones(len(self.range)) if i == 0 else np.concatenate((y, yi*np.ones(L)))
        y = np.tile(y,H)

        for i,zi in enumerate(self.depth):
            z = zi*np.ones(L*L) if i == 0 else np.concatenate((z, zi*np.ones(L*L)))
        
        return pd.DataFrame({"x" : x , "y" :  y , "z" : z,
                              "t" : t , "d" : d, "p" : p,
                              "u" : u , "v" : v, "w" : w})

