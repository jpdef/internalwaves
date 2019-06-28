#IW_FIELD
#Desc : Library for computing internal wave fields
#Auth : J. DeFilippis
#Date : 6-26-2019

import matplotlib.pyplot as plt

class iw_vmodes:
    """
    Desc: 
    A class to generate internal wave field modes via different
    numerical methods

    Attributes:
       depth : array
         A depth coordinate of choice meters,pressure etc
       range : array
         A depth coordinate of choice meters,pressure etc
       strat : func
         A function of stratification whose arguement is the 
         depth coordinate strat(depth)
    """

    def __init__(self,depth,N):
        pass

    
    def omegatok(omega,m):
        """
        Desc:
        Transforms frequency to horizontal wave number through 
        the dispersion relation
        Parameters:
            omega : frequency of the wave
            m     : the vertical mode number
        """
        pass 
