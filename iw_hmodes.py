#FFTE (fast fourier transform easy)
#Desc : Library for fast fourier transforms 
#Auth : J. DeFilippis
#Date : 4-03-2019

##Interpolations cause errors when transforming back and forth...

import numpy as np
from scipy import interpolate


class ftte:
    """
    Desc: 
    Class builds off of numpy's fast fourier transform library, but makes
    passing in functions easier by interpolating those functions for the
    user

    Attributes:
       N : int
         number of grid points to evaluate the function
       L : float
         domain of the function in real space       
    """

    def __init__(self,N=2,L=1):
        """
        Desc : iniates class with a set of grid points
        """
        self.N = N
        self.grid  = np.linspace(0,L,N)
        self.kgrid = np.linspace(0,2*np.pi,N)
    
    
    def transform(self,f_in_space):
        """
        Desc : fourier transforms a function and returns another
               function defined on the domain 2 pi
        Params:
               @f : is the function that is fourier transformed.
        """
        f_on_grid  = f_in_space(self.grid)
        f_on_kgrid   = np.fft.fft(f_on_grid)
        f_in_kspace = interpolate.interp1d(self.kgrid,f_on_grid)
        return f_in_kspace
    
    
    def inv_transform(self,f_in_kspace):
        """
        Desc : inverse fourier transforms a function and returns another
               function defined on the domain L
        Params:
               @f_in_kspace :  is the function that is inverse fourier 
                            transformed.
        """
        f_on_kgrid  = f_in_kspace(self.kgrid)
        f_on_grid   = np.fft.ifft(f_on_kgrid)
        f_in_space = interpolate.interp1d(self.grid,f_on_grid)
        return f_in_space
    
