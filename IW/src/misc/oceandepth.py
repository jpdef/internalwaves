#Desc: Defines a object that can hold all parameters of an vertical ocean profile, evenly spaced between  points.
#Auth : J. DeFilippis
#date : 2/22/19

import numpy as np
import scipy as sp
import seawater as sw

class odepth:
    """
    Desc:
        Depth object defines the vertical coordinate in multiple dimensions.It tracks 
        paramters such as pressure salinity and temperature at different depths. It assumes
        all paramters that are feed into the object were measured at equal distances.
    Fields:
        bottom : defines the maximum depth 
        res : defines the amount of discrete points that make the vertical coordinate
    """


    def __init__(self,max_depth,res,lat=15):
        self.max_depth = max_depth
        self.res  = res
        self.lat  = lat
        self.dims = dict()
        self.z = np.linspace(0,self.max_depth, self.res)
        self.dims['z'] = self.z


    def add_dim(self,dim_name,dim):
        """
        Desc :
            Interpolates a dimension, generates an array of equal size to z 
            and stores in dictionary
        """

        ztemp = np.linspace(0,self.max_depth,len(dim))
        f=sp.interpolate.interp1d(ztemp,dim)
        self.dims[dim_name] = f(self.z)


    def get_dim(self,dim_name):
        """
        Desc:
            Returns dimension from object
        """
        if dim_name not in self.dims:
            self.gen_dim(dim_name)
            return self.get_dim(dim_name)
        else:
            return self.dims[dim_name]

    def gen_dim(self,dim_name):
        """
        Desc:
            Auto-generates and stores other dimensions if they
            can be calculated from other stored dimensions
            ie:// density(S,T,P)

        """
        if dim_name is 'pressure' :
            z = self.get_dim('z')
            p = sw.eos80.pres(z,self.lat)
            self.add_dim('pressure',p)

        elif dim_name is 'density':
            s = self.get_dim('salinity')
            t = self.get_dim('temperature')
            p = self.get_dim('pressure')
            d = sw.eos80.dens(s,t,p)
            self.add_dim('density',d)

        elif dim_name is 'bfrq':
            s = self.get_dim('salinity')
            t = self.get_dim('temperature')
            p = self.get_dim('pressure')
            strat = sw.geostrophic.bfrq(s,t,p,lat=self.lat)
            self.add_dim('bfrq',np.hstack(strat[0]))  # Brunt Vaisala Frequency
            self.add_dim('ppv' ,np.hstack(strat[1]))  # Planetary Potential Vorticity
            self.add_dim('pmp' ,np.hstack(strat[2]))  # Midpoint pressure
        
        elif dim_name is 'svel':
            s = self.get_dim('salinity')
            t = self.get_dim('temperature')
            p = self.get_dim('pressure')
            svel = sw.svel(s,t,p)
            self.add_dim('svel',svel)

        else:
            print("Cannot be generated from stored information")
            return None
