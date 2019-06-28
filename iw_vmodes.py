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

    def __init__(self,depth,N):
        """
        Parameters:
          depth : array
              vertical coordinate of modes
          N : func, or array
              stratification of medium
        """
        self.depth = depth
        self.N = N
        self.modes = []
    

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
            self.modes = [ np.sin(np.pi*j*eta(self.depth) for j in range(0,len(self.depth) ]
            return np.sin(np.pi*j*eta(self.depth))
        
        elif isinstance (self.N, np.ndarray):
            No  = trapz(self.N,self.depth)
            eta = np.vectorize( lambda zn : (1/No)*trapz(self.N[zn:],self.depth[zn:]) )
            zn = np.arange(0,len(self.depth))
            self.modes = [ np.sin(np.pi*j*eta(zn) for j in range(0,len(self.depth) ]
            return np.sin(np.pi*j*eta(zn))
        
        else:
            print("N needs to be an array or function")
            return


    def gen_vmodes_evp(self,k=2*np.pi):
        """
        Desc : Solves helmoltz equation using an EVP solver 
               techique  
        Args :
               k : float 
                The horizontal wave number , default 2pi
        Returns:
               solution : array
                 A solution as a function of the depth coordinate
        """
        #Physical Properties
        H  = len(self.depth)
        delta = (self.depth[H-1] - self.depth[0])/H
        ksq = (k)**2
        
        #FDM stencil for centered second derivative
        D2  = self.gen_tri_fdm(H,delta)
        K2   = np.diag(np.ones(len(self.N))*ksq)
        N2   = np.diag(self.N)
        
        #Eigen value problem IW equation
        w,vr = eig(-N2,(D2-K2))
 
        self.modes = [ vr[:,m] for m in np.arange(0,len(vr)) ]
     
        return w,vr


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
