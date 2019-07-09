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

    def __init__(self,bfrq,freqs,alpha=0.5,num_modes=1,npoints=(100,50)):
        self.depth = np.linspace(0,alpha,npoints[0])
        self.range = np.linspace(0,1,npoints[1])
        self.vmodes = iwvm.iw_vmodes(bfrq,self.depth)
        self.num_modes = num_modes
        self.hwavenumbers = np.zeros([self.num_modes,len(freqs)])



    def construct_field(self):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        pass


    def construct_hwavenumbers(self,freqs):
        #Construct a set of dummy wave numbers from range
        hwavenumbers = [20*np.pi/(wl+0.1) for wl in self.range]
        
        #Make a interploations of dispersion for each mode
        self.construct_dispersion_curves(hwavenumbers)

        #Invert interpolation to compute wavenumber from input frequency
        self.hwavenumbers = self.transform_frequencies_to_wavenumber(freqs) 


    def construct_dispersion_curves(self,hwavenumbers):
        """
        Desc:
        Interpolates a dispersion curve from the eigenvalues of
        the vertical mode solver. 
        """
        self.dispcurves=[]
        for m in range(self.num_modes):
            freqs = []
            for wn in hwavenumbers:
                cp,vr = self.vmodes.gen_vmodes_evp(wn)
                freqs.append( np.sqrt(cp[m].real)*wn )
            f = interpolate.interp1d(hwavenumbers,freqs)
            self.dispcurves.append(f)


    def transform_frequencies_to_wavenumber(self,freqs):
        """
        Desc:
        Transforms frequency to horizontal wave number through 
        the dispersion relation
        Parameters:
            omega : frequency of the wave
            m     : the vertical mode number
        """
        hwavenumbers = np.zeros([self.num_modes,len(freqs)])
        for m in range(self.num_modes):
            for f,freq in enumerate(freqs):
                hwavenumbers[m,f] = self.invert_for_frequency(freq,m)
        return hwavenumbers


    def invert_for_frequency(self,freq,m):
        """
        Desc:
        Transforms an individual frequncy to a wavenumber via
        dispersion curve that has been generated
        """
        krange = [20*np.pi/(wl+0.1) for wl in self.range]
        if self.dispcurves is None:
           self.construct_dispersion_curves()
        f =  self.dispcurves[m]
        fshift = lambda x : f(x) - freq
        try:
            sol = scipy.optimize.fsolve(fshift,krange[-1])
            print(sol)
            return sol
        except ValueError:
            return -1
        
