#Desc : Library for constructing linear models
#Auth : J. DeFilippis
#Date : 7-25-2020

import numpy as np 
from ModeSolver import InternalWaveModes
from ModeSolver import QuasiGeostrophicModes
from tqdm import tqdm
from functools import partial

class LinearModel:
    """
    Desc:
    A class to make general linear models with matrices
    @ basis (function) takes parameters produce a second function
            that is evaluate at a coordinate
    @ parameters (array) variables for the basis function e.g. k,l,m
    """
    
    def __init__(self, basis , parameters):
        """
        Takes a vector of parameters to form columns of matrix
        and a vector of coordinates to form rows using the basis
        generator function
        """
        self.basis         = basis
        self.parameters    = parameters
        self.columns = sum([len(x) for x in parameters])
        
    def __add__(self,other):
        #Concatenate columns
        sumlm = LinearModel(self.basis+other.basis,self.parameters + other.parameters)     
        sumlm.columns = self.columns + other.columns
        return sumlm
    
    
    def model_matrix(self,coordinates):
        rows = len(coordinates)
        H = np.zeros(shape=(rows,self.columns))
        cn = 0 
        for i, b in enumerate(self.basis):
            for j,p in enumerate(self.parameters[i]):
                f = b(*p)
                H[:,cn] =  f(coordinates) 
                cn +=1 
        
        return H

    
class InternalWaveModel(LinearModel):
    def __init__(self,parameters,N2,Z,scalar_gradient=None,free_surface=False,
                      latitude=30):
        self.frequencies = np.unique(parameters['f'])
        
        
        self.iwmodes = self.generate_modes(N2,Z,scalar_gradient,free_surface,
                                           latitude)
        bases = [partial(iw_modal_sine,iw_modes=self.iwmodes),
                 partial(iw_modal_cosine,iw_modes=self.iwmodes)]
        
        super().__init__(bases,[parameters]*2)
    
    def generate_modes(self,N2,Z,scalar_gradient,free_surface,latitude):
        iwmodes = []
        for frequency in self.frequencies:
            m = InternalWaveModes(N2,Z,frequency,scalar_gradient=scalar_gradient,
                                  free_surface=free_surface,
                                  latitude=latitude)
            iwmodes.append(m)
        return iwmodes
    
    
class QuasiGeostrophicModel(LinearModel):
    def __init__(self,parameters,N2,Z,scalar_gradient=None,free_surface=False,
                      latitude=30):
        self.frequencies = np.unique(parameters['f'])
        self.qgmodes = self.generate_modes(N2,Z,scalar_gradient,free_surface,
                                           latitude)
        bases = [partial(iw_modal_sine,iw_modes=self.qgmodes),
                 partial(iw_modal_cosine,iw_modes=self.qgmodes)]
        
        super().__init__(bases,[parameters]*2)
    
    def generate_modes(self,N2,Z,scalar_gradient,free_surface,latitude):
        qgmodes = []
        for frequency in self.frequencies:
            m = QuasiGeostrophicModes(N2,Z,frequency,scalar_gradient=scalar_gradient,
                                  free_surface=free_surface,
                                  latitude=latitude)
            qgmodes.append(m)
        return qgmodes

        
        
    
class BasisFunctions:

    def __init__(self, function=None,**kwargs):
        self.kwargs = kwargs
        if function:
            self.func = getattr(self,function)(kwargs)
        
    def __mul__(self,other):
        l = len(self.kwargs)
        def function(*args):
            return self.func(*args[:l])*other.func(*args[l:])
        b = BasisFunctions()
        b.kwargs = {**self.kwargs, **other.kwargs}
        b.func = function
        return b

    
    def test(self, kwargs):
        print(len(kwargs) )
        
        
    def cosine(self, kwargs):
        m = kwargs['m']
        def function(*args):
            x = args[0]
            print('func2: ',np.cos(2*np.pi*m*x))
            return  np.cos(2*np.pi*m*x)
        
        return function
    
    def cosine_2d(self, kwargs):
        k = kwargs['k']
        l = kwargs['l']
        def function(*args):
            x = args[0];y=args[1];
            arg = 2*np.pi*(k*x + l*y)
            return np.cos(arg)
        
        return function
    
    def cosine_2d_time(self, kwargs):
        k = kwargs['l']
        l = kwargs['k']
        o = kwargs['o']
        def function(*args):
            x = args[0];y=args[1];t=args[2]
            arg = 2*np.pi*(k*x + l*y - o*t)
            print('func1:',np.cos(arg))
            return np.cos(arg)
        
        return function
    
    
def iw_modal_cosine(o,m,t,iw_modes):
    #print('modal_cosine')
    freqs = np.array( [mode.frequency for mode in iw_modes] )
    i = np.argmin(abs( o - freqs ) )
    K =  iw_modes[i].hwavenumber(m).real
    k =  K*np.cos(t*np.pi/180)
    l =  K*np.sin(t*np.pi/180)
    
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        ret = iw_modes[i].evaluate_mode(m,args['z'])*np.cos(arg)
        return ret.real
    
    return function


def iw_modal_sine(o,m,t,iw_modes):
    freqs = np.array( [mode.frequency for mode in iw_modes] )
    i = np.argmin(abs( o - freqs ) )
    K =  iw_modes[i].hwavenumber(m).real
    k =  K*np.cos(t*np.pi/180)
    l =  K*np.sin(t*np.pi/180)
    
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        ret = iw_modes[i].evaluate_mode(m,args['z'])*np.sin(arg)
        return ret.real
    
    return function

def modal_plane_wave_c(K,t,m,o,modes):
    mode = modes[m]
    k =  K*np.cos(t*np.pi/180)
    l =  K*np.sin(t*np.pi/180)
    
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        return mode(args['z'])*np.cos(arg)
    
    return function
      
    
    
def modal_plane_wave_s(K,t,m,o,modes):
    mode = modes[m]
    k =  K*np.cos(t*np.pi/180)
    l =  K*np.sin(t*np.pi/180)
    
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        return mode(args['z'])*np.sin(arg)
    
    return function


def modal_cosine(k,l,m,o,modes):
    mode = modes[m]
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        return mode(args['z'])*np.cos(arg)
    
    return function
      
    
    
def modal_sine(k,l,m,o,modes):
    mode = modes[m]
    
    def function(args):
        arg = 2*np.pi*(k*args['x'] + l*args['y']- o*args['t'])
        return mode(args['z'])*np.sin(arg)
    
    return function

