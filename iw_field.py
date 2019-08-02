#IW_FIELD
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
                 freqs=np.array([]), hwavenumbers=np.array([]),
                 weights=np.array([]), bfrq=np.array([]), 
                 randomphase=True,offset=0,
                 modes=np.array([1]),npoints=(100,50),
                 scalarfield=np.array([])):
        
        if not freqs.size and not hwavenumbers.size:
            print("""Internal Wave field needs a set of frequencies 
                   or horizontal wavenumbers""")
            return
        else:
            print("Intializing wavefield")
        
        #Set all fields in object
        self.set_attributes(bfrq,iwrange,iwdepth,modes,npoints,randomphase,offset)        
        
        #Compute phase speeds and vertical structure functions
        self.init_dispersion(freqs,hwavenumbers,weights)

        #Compute 2D field from phase speeds and structure functions
        self.field = self.construct_field(freqs,hwavenumbers)
        
        #Impose IWF on scalarfield (e.g. temperature)
        if scalarfield.size:
            self.scalarfield = scalarfield
            self.scalarfield.distort(self.field)
        else:
            self.scalarfield = scalarfield

    
    def construct_field(self,freqs,hwavenumbers):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        if freqs.size and not hwavenumbers.size:
           return self.construct_field_from_freqs()
        elif hwavenumbers.size and not freqs.size:
           return self.construct_field_from_wavenumbers()

    
    def construct_field_from_freqs(self,new_weights=np.array([])):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        field = np.zeros(shape=(len(self.depth),len(self.range)),dtype=complex)
        if not self.field_components:
            self.construct_field_components(self.freqs.shape[1]) 
        
        for n,f in enumerate(self.field_components):
           field += new_weights[n]*f if new_weights.size else self.weights[n]*f
        
        return field   

 
    def construct_field_from_wavenumbers(self):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        field = np.zeros(shape=(len(self.depth),len(self.range)),dtype=complex)
        if not self.field_components:
            self.construct_field_components(self.hwavenumber.shape[1]) 

        for n,f in enumerate(self.field_components):
           field += self.weights[n]*f
 
        return field

  
    def construct_field_components(self,number_of_components):
        """
        Desc:
        Fills list of 2D array with the modal components that sum to 
        make the total field. Components are stored for efficiency
        """
        for n in range(number_of_components):
            self.field_components.append(self.construct_field_component(n))


    def construct_field_component(self,n):
        """
        Desc:
        Constructs a 2D wave field componenent for a specific frequency/wavenumber
        depending on the choosen dependent variable. These components are then summed
        to give the total wavefield.
        """
        zeta = np.zeros(shape=(len(self.depth),len(self.range)),dtype=complex)
        for i,m in enumerate(self.modes):
            k     = self.hwavenumbers[i,n]
            phi   = self.vertical_comp[:,m,n]
            if self.randomphase:
                rp  = 2*np.pi*np.random.rand()
                psi = np.exp(2*np.pi*1j * k * (self.range - self.offset) + rp ) 
            else:
                psi = np.exp(2*np.pi*1j * k * (self.range - self.offset)  )
            zeta += np.outer(phi,psi)
        return ( zeta / len(self.modes) ) 


    def init_dispersion(self,freqs,hwavenumbers,weights):
        """
        self.vertical_component[z_n, m, f_n] 
            z_n depth grid point for vertical functions
            m   mode number
            f_n associated frequency
        """
        D = len(self.depth)
        if freqs.size and not hwavenumbers.size:
            print("Field with input of frequencies")
            nf = len(freqs)
            self.vertical_comp = np.ndarray(shape=(D,D,nf ) )
            self.freqs = np.tile(freqs,( len(self.modes) , 1))
            self.hwavenumbers = self.construct_hwavenumbers(freqs)
            self.weights = weights if weights.size else np.ones(nf)
        elif hwavenumbers.size and not freqs.size:
            print("Field with input of wavenumbers")
            nk = len(hwavenumbers)
            self.vertical_comp = ndarray(shape=(D,D,nk) )
            self.hwavenumbers = np.tile(hwavenumbers,( len(self.modes) ,1))
            self.freqs = self.construct_frequencies(hwavenumbers)
            self.weights = weights if weights.size else np.ones(nk)


    def construct_hwavenumbers(self,freqs):
        """
        Desc:
        Constructs a set of horizontal wavenumbers from the
        set of input frequencies via the dispersion relations
        """
        #Make a interploations of dispersion for each mode
        self.construct_dispersion_curves(True,freqs)

        #Invert interpolation to compute wavenumber from input frequency
        return self.transform_frequencies_to_wavenumber(freqs) 
    

    def construct_frequencies(self,hwavenumbers):
        """
        Desc:
        Constructs a set of frequencies from the set of input frequencies
        input horizontal wavenumbers  via the dispersion relations
        """
        #Make a interploations of dispersion for each mode
        self.construct_dispersion_curves(False,hwavenumbers)
        
        #Invert interpolation to compute wavenumber from input frequency
        return self.transform_wavenumber_to_frequencies(hwavenumbers) 


    def construct_dispersion_curves(self,isfreqs,independent_vars):
        """
        Desc:
        Interpolates a dispersion curve from the eigenvalues of
        the vertical mode solver. 
        """
        self.dispcurves=[]
        for m in self.modes:
            dependent_vars = []
            for itr,ivar in enumerate(independent_vars):
                if isfreqs:
                    cp,vr = self.vmodes.gen_vmodes_evp_f(ivar)
                    dependent_vars.append( ivar/np.sqrt( cp[m]))
                    self.vertical_comp[:,:,itr] = vr
                else:
                    cp,vr = self.vmodes.gen_vmodes_evp_k(ivar)
                    dependent_vars.append( np.sqrt(cp[m])*ivar )
                    self.vertical_comp[:,:,itr] = vr
                    
            #f = interpolate.interp1d(independent_vars,dependent_vars,kind='cubic')
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

    def update_field(self,new_weights):
        """
        Desc:
        Efficient way to timestep or change distrubution after having calculated all
        the horizontal and vertical components
        """
        self.field = self.construct_field_from_freqs(new_weights)
        if self.scalarfield.size:
            self.scalarfield.distort(self.field)

 
    def set_attributes(self,bfrq,iwrange,iwdepth,modes,npoints,
                       randomphase,offset):
        """
        Desc:
        Helper function for constructor to set all the various fields
        """
        self.range = iwrange
        self.depth = iwdepth 
        self.bfrq = bfrq if bfrq.size else self.cannonical_bfrq()
        self.vmodes = iwvm.iw_vmodes(self.depth,self.bfrq)
        self.modes = modes
        self.randomphase = randomphase
        self.offset = offset
        self.field_components = [] 
        self.field = np.zeros(shape=(len(self.depth),len(self.range)),dtype=complex)


    def cannonical_bfrq(self):
        """
        Desc:
        Cannonical example of a mid ocean Brunt Viasala frequency
        depth profile. Mainly for testing purposes
        """
        d = max(self.depth)
        sigma = 22 + 2.5*np.tanh(2*np.pi*(self.depth-.15*d)/d)
        N     = np.sqrt(np.gradient(sigma))/5.0
        cph   = 3600/(2*np.pi)
        return N/cph


    def to_dataframe(self):
        """
        Desc:
        Converts the 2D array internal wave field to a panda's
        dataframe for use and convience
        """
        zeta = self.field.real.flatten()
        x    = np.ndarray(shape=(1,),dtype=float)
        z    = np.ndarray(shape=(1,),dtype=float)
       
        for i,zi in enumerate(self.depth):
            x = self.range if i==0 else np.concatenate( (x, self.range) )
           
        for i,zi in enumerate(self.depth):
            if (i==0): 
                z = zi*np.ones(len(self.range)) 
            else: 
                z = np.concatenate( (z, zi*np.ones(len(self.range)) ) )
        
        return pd.DataFrame({"range" : x , "depth" : z , "disp" : zeta})
