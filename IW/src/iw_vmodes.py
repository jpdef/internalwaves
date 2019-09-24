#IW_VMODES
#Desc : Library for solving vertical modes of internal waves
#Auth : J. DeFilippis
#Date : 2-15-2019

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp, quad,trapz
from scipy.linalg import eig ,inv

class iw_vmodes:
    """
    Desc: 
    A class to generate internal wave vertical modes via different
    numerical methods

    Attributes:
       depth : array
         A depth coordinate of choice meters,pressure etc
       strat : func
         A function of stratification whose arguement is the 
         depth coordinate strat(depth)
    """

    def __init__(self,depth,N=np.array([]),T=np.array([]),f=0):
        """
        Parameters:
          depth : array
              vertical coordinate of modes
          N : func, or array
              stratification of medium
        """
        #Set Attributes
        self.set_attributes(depth,N,T,f)

        #Generate Modes 
        lamb,vr = self.gen_vmodes_evp()
        self.modes        = [ vr[:,m] for m in np.arange(0,len(vr)) ]
        self.hwavenumbers =  f / np.sqrt(lamb)

        #Generate other parameters
        self.int_modes = self.gen_integral_modes()


    def set_attributes(self,depth,N=np.array([]),T=np.array([]),f=0):
        """
        Desc : Set various attributes for the class
        """
        self.depth = depth
        if N.size: 
            self.N = N  
        else if T.size :
            self.N = self.compute_bfrq(T,S)
        else :
            self.N = self.cannonical_bfrq()
        
        self.freq = f if f > 0 else 2*np.pi/(3600*24)
        self.modes = []
        self.hwavenumbers = []
    

    def gen_vmodes_wkb(self,j):
        """
        Desc : Solves helmholtz equation using WKB
               approximation
        Params :
          j : int
              mode number
        Returns:
               solution : array
                 A solution as a function of the depth coordinate

        """
        
        #Integrate over stratification
        if callable(self.N):
            No  = quad(self.N,self.depth[-1],self.depth[0])[0]
            eta = np.vectorize(lambda z :  (1/No)*quad(self.N,self.depth[-1],z)[0])
            self.modes = [ np.sin(np.pi*j*eta(self.depth)) for j in np.arange(0,len(self.depth)) ]
            return np.sin(np.pi*j*eta(self.depth))
        
        elif isinstance (self.N, np.ndarray):
            No  = trapz(self.N,self.depth)
            eta = np.vectorize( lambda zn : (1/No)*trapz(self.N[zn:],self.depth[zn:]) )
            zn = np.arange(0,len(self.depth))
            self.modes = [ np.sin(np.pi*j*eta(zn)) for j in np.arange(0,len(self.depth)) ]
            return np.sin(np.pi*j*eta(zn))
        
        else:
            print("N needs to be an array or function")
            return


 
    def gen_vmodes_evp(self,f=2*np.pi):
        """
        Desc : Solves helmoltz equation using an EVP solver 
               techique  
        Args :
               f : float 
                 Radial frequency , default 2pi
        Returns:
               solution : array
                 A solution as a function of the depth coordinate
        """
        #Physical Properties
        H  = len(self.depth)
        delta = (self.depth[H-1] - self.depth[0])/H
        fsq = (f)**2
        
        #FDM stencil for centered second derivative
        D2  = self.gen_tri_fdm(H,delta)
        F2  = np.diag(np.ones(len(self.N))*fsq)
        N2  = np.diag(self.N)
        
        #Eigen value problem IW equation
        lamb,vr = eig((F2-N2),D2)
        
        return lamb,vr
    
        
    def gen_integral_modes():
        """
        Desc : 
        Generates the integral of each mode so that variables such as
        pressure can be computed
        """
        int_modes = []
        for m in self.modes:
            chi = 1025 * ((self.N**2 - omega**2)/omega) * m  
            int_mode  = trapz(self.chi[0:i],self.depth[0:i] for i in len(chi))
            int_modes.append(int_mode)
        
        return int_modes


    #Normalization might now work for model
    def normalize(self,vr):
        """
        Desc : 
        Normalizes each mode structure by taking the potential energy
        depth integral to be 1
        """
        for m in range(vr.shape[1]):
            A = np.sqrt( 1.0/np.sum((vr[:,m]**2)*(self.N)**2) )
            vr[:,m] = A*vr[:,m]
        return vr

    
    #Consider making this more generic
    def gen_tri_fdm(self,N,delta):
        """
        Desc : Generates a tri-diagonal matrix for 2nd order 
               derivative using a finite difference
        Returns :
                matrix : np matrix
        """
        ret = np.zeros(shape=(N,N))
        for i in range(1,N-1):
            ret[i,i]   = -2 
            ret[i,i-1] = 1 
            ret[i,i+1] = 1 
        
        ret[0,0]   = ret[N-1,N-1] = -2
        ret[0,1]   = ret[N-1,N-2] = 1

        ret /= delta**2

        return ret

    
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


    def mode_plot(self,n=3):
        """
        Desc : Generates a depth plot of stratriciation 
               against the first n modes
        Returns :
                fig,ax : matplotlib figure and axes
        """
        f, (ax1,ax2) = plt.subplots(1,2,sharey=True)
        
        #Stratification
        ax1.invert_yaxis()
        ax1.plot(3600*np.sqrt(self.N/(2*np.pi)),self.depth  )
        ax1.set_title("Stratification")
        ax1.set_xlabel("CPH")
    
        #Modes
        for i in np.arange(0,n):
            m = self.modes[i]
            ax2.plot(m, self.depth, label="mode "+str(m))
        ax2.set_title("Modes")
       
        return f, (ax1,ax2) 
