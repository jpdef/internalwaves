#Desc : Library for computing internal wave fields
#Auth : J. DeFilippis
#Date : 6-26-2019

import matplotlib.pyplot as plt
import iw_vmodes as iwvm
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
                 offset=[0,0]):
        
        print("Intializing wavefield")
        
        #Set all fields in object
        self.set_attributes(bfrq,iwrange,iwdepth,modes,offset)        
        
        #Compute phase speeds and vertical structure functions
        self.init_dispersion(freqs,amplitudes)

        #Compute 3D field from phase speeds and structure functions
        self.disp_field = self.construct_disp_field(freqs)
        self.velo_field = self.construct_velo_field(freqs)
        
    
    def construct_disp_field(self,step=np.array([])):
        """
        Desc:
        Constructs a 3D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        field = self.empty_field()
                
        #Initialize Field Components 
        if not self.disp_field_components:
            for n in range(self.nfreqs):
                self.disp_field_components.append(self.construct_field_component(n))
        
        #Update Field Component Values
        for n,fc in enumerate(self.disp_field_components):
           field += step[n]*fc if step.size else fc 
        
        return field   
   
    #TODO: need to multiply k/|k| and l/|k| for each horiztonal component 
    def construct_velo_field(self,step=np.array([])):
        """
        Desc:
        Constructs a 3D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        field = self.empty_field()
                
        #Initialize Field Components 
        if not self.velo_field_components:
            for n in range(self.nfreqs):
                self.velo_field_components.append(self.construct_velo_field_component(n))
        
        #Update Field Component Values
        for n,fc in enumerate(self.velo_field_components):
           field += step[n]*fc if step.size else fc 
        
        return field   

 
    def construct_disp_field_component(self,n):
        """
        Desc:
        Constructs a 3D wave field componenent for a specific frequency/wavenumber
        depending on the choosen dependent variable. These components are then summed
        to give the total wavefield.
        """
        zeta = self.empty_field()
       
        for i,m in enumerate(self.modes):
            psi  = self.horizontal_comp(i,n)
            phi  = self.vertical_comp[:,m,n]
            for zn in range(len(self.depth)):
                zeta[zn,:,:] += psi*phi[zn]
       
        return zeta 
   
    
    def construct_velo_field_component(self,n):
        """
        Desc:
        Constructs a 3D wave field componenent for a specific frequency/wavenumber
        depending on the choosen dependent variable. These components are then summed
        to give the total wavefield.
        """
        zeta = self.empty_field()
       
        for i,m in enumerate(self.modes):
            psi    = self.horizontal_comp(i,n)
            phiz   = self.vertical_grad[:,m,n]
            for zn in range(len(self.depth)):
                zeta[zn,:,:] += psi*phiz[zn]
       
        return zeta 
   
    
    def horizontal_comp(self,i,n):
        """
        Desc:
        Constructs a plane wave in the horizontal with a specific 
        wavenumber and heading  
        """
        xx,yy = np.meshgrid(self.range - self.offset[0],self.range - self.offset[1])
        kmag  = self.hwavenumbers[i,n]
        amps  = self.get_amplitude(i,n) 
        psi   = np.zeros ( shape=xx.shape, dtype='complex128' )
        for ah in list(zip( amps['amps'],amps['headings'])):
            psi  += self.plane_wave(kmag,ah[0],ah[1],xx,yy) 
        return(psi)
   
    
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
        field = np.zeros(shape=(len(self.range),len(self.range),len(self.depth)),dtype=complex)
        return field


    def default_amplitudes(self,n):
        return([{'ar': 1 ,'ai' : 0 ,'theta' : np.pi/4} for n in range(n)])


    def init_dispersion(self,freqs,amplitudes):
        """
        self.vertical_comp[z_n, m, f_n] 
        self.vertical_grad[z_n, m, f_n] 
            z_n depth grid point for vertical functions
            m   mode number
            f_n associated frequency
        """
        D = len(self.depth)
        self.nwaves = len(freqs)*len(self.modes)
        self.nfreqs = len(freqs)
        
        self.vertical_comp = np.ndarray(shape=(D,D,self.nwaves ) )
        self.vertical_grad = np.ndarray(shape=(D,D,self.nwaves ) )
        self.freqs         = np.tile(freqs,( len(self.modes) , 1))
        self.hwavenumbers  = self.construct_hwavenumbers(freqs)
        self.amplitudes    = amplitudes if amplitudes else self.default_amplitudes(self.nwaves)


    def construct_hwavenumbers(self,freqs):
        """
        Desc:
        Constructs a set of horizontal wavenumbers from the
        set of input frequencies via the dispersion relations
        """
        #Make a interploations of dispersion for each mode
        self.construct_dispersion_curves(freqs)

        #Invert interpolation to compute wavenumber from input frequency
        return self.transform_frequencies_to_wavenumber(freqs) 
    

    def construct_dispersion_curves(self,independent_vars):
        """
        Desc:
        Interpolates a dispersion curve from the eigenvalues of
        the vertical mode solver. In this case independent_var
        is the frequencies and dependent_var is wavenumber 
        """
        self.dispcurves=[]
        for m in self.modes:
            dependent_vars = []
            for itr,ivar in enumerate(independent_vars):
                cp,vr = self.vmodes.gen_vmodes_evp_f(ivar)
                dependent_vars.append( ivar/np.sqrt( cp[m]))
                self.vertical_comp[:,:,itr] = vr
                self.vertical_grad[:,m,itr] = np.gradient(vr[:,m],self.zdiff)
            
            #Using dictionary instead of interpolation so we can have small sets
            #of frequencies that the interpolation cannot support
            d  = dict(zip (independent_vars, dependent_vars ))
            self.dispcurves.append(d)
   
 
    def transform_wavenumber_to_frequencies(self,hwavenumbers):
        """
        Desc:
        Converts  horizontal wavenumbers to frequencies using
        the precomputed dispersion curves
        """
        freqs = []
        if self.dispcurves:
            for d in self.dispcurves:
                freqs.append([ d[h] for h in hwavenumbers ])
            return np.vstack(freqs)
        else:
           raise KeyError
   
 
    def transform_frequencies_to_wavenumber(self,freqs):
        """
        Desc:
        Converts frequencies to horizontal wavenumbers using
        the precomputed dispersion curves
        """
        hwavenumbers = []
        if self.dispcurves:
           for d in self.dispcurves:
               hwavenumbers.append( [ d[f] for f in freqs] )    
           return np.vstack(hwavenumbers)
        else:
           raise KeyError


    def update_field(self,step):
        """
        Desc:
        Efficient way to timestep or change distrubution after having calculated all
        the horizontal and vertical components
        """
        self.disp_field = self.construct_disp_field(step)

 
    def set_attributes(self,bfrq,iwrange,iwdepth,modes,offset):
        """
        Desc:
        Helper function for constructor to set all the various fields
        """
        self.range   = iwrange
        self.depth   = iwdepth
        self.zdiff   = np.average(np.diff(iwdepth)) 
        self.bfrq    = bfrq if bfrq.size else self.cannonical_bfrq()
        self.vmodes  = iwvm.iw_vmodes(self.depth,self.bfrq)
        self.modes   = modes
        self.offset  = offset
        self.disp_field_components = [] 
        self.disp_field = np.zeros(shape=(len(self.range),len(self.range),len(self.depth)),
                              dtype=complex)

    #Add the depth spacing to the stratifcation grad and check if this works for R script
    def cannonical_bfrq(self):
        """
        Desc:
        Cannonical example of a mid ocean Brunt Viasala frequency
        depth profile. Mainly for testing purposes
        """
        d = max(self.depth)
        sigma = 22 + 2.5*np.tanh(2*np.pi*(self.depth-.15*d)/d)
        N     = np.sqrt(np.gradient(sigma))
        return N/3600 #convert to cycle per second


    def to_dataframe(self,coords=[],time=0):
        """
        Desc:
        Converts the 3D array internal wave field to a panda's
        dataframe for use and convience
        """
    
        return self.select_data(coords,time) if coords else self.flatten_data() 
    
    def select_data(self,coords,time):
         x = [self.range[c[0]] for c in coords] 
         y = [self.range[c[1]] for c in coords] 
         z = [self.depth[c[2]] for c in coords]  
         Z = [self.disp_field.real[c[2],c[1],c[0]] for c in coords]
         t = np.repeat(time,len(x))
        
         return pd.DataFrame({"x" : x , "y" :  y , "z" : z , "t" : t, "disp" : Z})
    

    def flatten_data(self):
        zeta = self.disp_field.real.flatten()
        x    = np.ndarray(shape=(1,),dtype=float)
        y    = np.ndarray(shape=(1,),dtype=float)
        z    = np.ndarray(shape=(1,),dtype=float)
       
        L = len(self.range)
        H = len(self.depth)

        x = np.tile(self.range,L*H) 
           
        for i,yi in enumerate(self.range):
            y = yi*np.ones(len(self.range)) if i == 0 else np.concatenate((y, yi*np.ones(L)))
        y = np.tile(y,H)

        for i,zi in enumerate(self.depth):
            z = zi*np.ones(L*L) if i == 0 else np.concatenate((z, zi*np.ones(L*L)))
        
        return pd.DataFrame({"x" : x , "y" :  y , "z" : z , "disp" : zeta})

