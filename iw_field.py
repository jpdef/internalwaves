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
        #Construct a set of dummy wave numbers from range
        self.hwavenumbers = [20*np.pi/(wl+0.1) for wl in self.range] 


    def construct_field(self):
        """
        Desc:
        Constructs a 2D wave field from the vertical and horizontal
        components by taking an outer product of the two vectors
        """
        pass


    def construct_dispersion_curves(self):
        """
        Desc:
        Interpolates a dispersion curve from the eigenvalues of
        the vertical mode solver. 
        """
        self.dispcurves=[]
        for m in range(self.num_modes):
            omegas = []
            for wn in self.hwavenumbers:
                cp,vr = self.vmodes.gen_vmodes_evp(wn)
                omegas.append( np.sqrt(cp[m].real)*wn )
            f = interpolate.interp1d(self.hwavenumbers,omegas)
            self.dispcurves.append(f)

    def transform_to_k(omega,m):
        """
        Desc:
        Transforms frequency to horizontal wave number through 
        the dispersion relation
        Parameters:
            omega : frequency of the wave
            m     : the vertical mode number
        """
        pass 
