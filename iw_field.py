#IW_FIELD
#Desc : Library for computing internal wave fields
#Auth : J. DeFilippis
#Date : 6-26-2019

import matplotlib.pyplot as plt
import iw_vmodes as iwvm
import scipy
from scipy import interpolate
import numpy as np

class InternalWaveField:
    """
    Desc: 
    A class to generate internal wave field modes via different
    numerical methods

    """

    def __init__(self,freqs=[],hwavenumbers=np.array([]),
                 bfrq=np.array([]), alpha=0.5,num_modes=1,npoints=(100,50)):
        
        if not freqs.size and not hwavenumbers.size:
            print("""Internal Wave field needs a set of frequencies 
                   or horizontal wavenumbers""")
            return
        else:
            print("Intializing wavefield")
        
        #Set all fields in object
        self.set_attributes(bfrq,alpha,num_modes,npoints)        
        
        #Compute phase speeds and vertical structure functions
        self.init_dispersion(freqs,hwavenumbers)

        #Compute 2D field from phase speeds and structure functions
        self.construct_field()
    

    def construct_field(self):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        pass

    
    def set_attributes(self,bfrq,alpha,num_modes,npoints):
        self.depth = np.linspace(0,alpha,npoints[0])
        self.range = np.linspace(0,1,npoints[1])
        self.bfrq = bfrq if bfrq.size else self.cannonical_bfrq()
        self.vmodes = iwvm.iw_vmodes(self.depth,self.bfrq)
        self.num_modes = num_modes


    def init_dispersion(self,freqs,hwavenumbers):
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
            self.freqs = np.tile(freqs,(self.num_modes,1))
            self.hwavenumbers = self.construct_hwavenumbers(freqs)
        elif hwavenumbers.size and not freqs.size:
            print("Field with input of wavenumbers")
            nk = len(hwavenumbers)
            self.vertical_comp = ndarray(shape=(D,D,nk) )
            self.hwavenumbers = np.tile(hwavenumbers,(self.num_modes,1))
            self.freqs = self.construct_frequencies(hwavenumbers)


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
        for m in range(self.num_modes):
            dependent_vars = []
            for itr,ivar in enumerate(independent_vars):
                if isfreqs:
                    cp,vr = self.vmodes.gen_vmodes_evp_f(ivar)
                    #dependent_vars.append( ivar/np.sqrt( abs(cp[m].real)))
                    dependent_vars.append( ivar/np.sqrt( cp[m]))
                    self.vertical_comp[:,:,itr] = vr
                else:
                    cp,vr = self.vmodes.gen_vmodes_evp_k(ivar)
                    dependent_vars.append( np.sqrt(cp[m].real)*ivar )
                    self.vertical_comp[:,:,itr] = vr
                    
            f = interpolate.interp1d(independent_vars,dependent_vars,kind='cubic')
            self.dispcurves.append(f)
     
 
    def transform_wavenumber_to_frequencies(self,hwavenumbers):
        if self.dispcurves:
            freqs = np.vstack( [ f(hwavenumbers) for f in self.dispcurves ])
            return freqs
        else:
           raise ValueError
   
 
    def transform_frequencies_to_wavenumber(self,freqs):
        if self.dispcurves:
           hwavenumbers = np.vstack( [ f(freqs) for f in self.dispcurves ])
           return hwavenumbers
        else:
           raise ValueError

  
    def cannonical_bfrq(self):
        """
        Desc:
        Cannonical example of a mid ocean Brunt Viasala frequency
        depth profile. Mainly for testing purposes
        """
        d = max(self.depth)
        sigma = 22 + 2.5*np.tanh(2*np.pi*(self.depth-.15*d)/d)
        N     = np.sqrt(np.gradient(sigma))/5.0
        return N         
