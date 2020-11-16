#Desc : Library for solving modes
#Auth : J. DeFilippis
#Date : 7-25-2020

import numpy as np
from scipy.linalg import eig

class ModeSolver:
    """
    One dimensional eigen value solver in  the form of 
    LH v = - e D2v 
    @ domain (array) that defines the points to solve on
    @ left_matrix (array) matrix on the left hand side of the eq.
    @ boundary (tuple) values at enpoints of domain
    @ free_surface (bool) free surface boundary at top v' = v
    """
    def __init__(self,domain,left_matrix=None,
                 boundary=[0,0],free_surface=False):
        
        
        #Check that domain is regular
        if (max(np.diff( np.diff(domain) ) ) > 1e-5 ):
            print("Domain is not regular by 1e-5")
            raise
            
        if (left_matrix is None):
            left_matrix = -1*np.identity(len(domain))
    
        #Check that right matrix has correct dimensions
        elif (left_matrix.shape[0] != left_matrix.shape[1]):
            print("Right hand matrix not-square")
            raise
            
        elif (left_matrix.shape[0] != len(domain)):
            print("Right hand matrix dim != domain dim")
            raise
        
        self.domain = domain
        self.npts = len(domain) - 2
        self.delta = np.mean(np.diff(domain))
        self.LHM = left_matrix
        self.bounds = boundary
        self.free_surface = free_surface
        self.mode_funcs = None
        
        
    def solve(self):
        if self.free_surface:
            D2 = self.centered_2nd_derivative_fs()
            LH = self.LHM[:-1,:-1]
            eigval,eigvec = eig(LH,-D2)
            eigvec = np.vstack([eigvec,
                               self.bounds[1]*np.ones(self.npts+1)])
            
        else:
            D2 = self.centered_2nd_derivative()
            LH = self.LHM[1:-1,1:-1]
            eigval,eigvec = eig(LH,-D2)
            eigvec = np.vstack([self.bounds[0]*np.ones(self.npts),
                               eigvec,
                               self.bounds[1]*np.ones(self.npts)])
            
        ind = np.argsort(eigval)[::-1]
        return [eigval[ind],eigvec[:,ind]]
    
    
    def centered_2nd_derivative_fs(self):
        """
        Desc : Generates a tri-diagonal matrix for 2nd order 
               derivative using a finite difference
        Returns :
                matrix : np matrix
        """
        N = self.npts + 1 
        D2 = np.zeros(shape=(N,N))
        for i in range(1,N-1):
            D2[i,i]   = -2 
            D2[i,i-1] = 1 
            D2[i,i+1] = 1 
        
        D2[0,1]  = -9.8*self.delta/2 
        D2[N-1,N-1] = -2
        D2[N-1,N-2] = 1
        
        D2 /= self.delta**2 
        return D2
        
        
    def centered_2nd_derivative(self):
        """
        Desc : Generates a tri-diagonal matrix for 2nd order 
               derivative using a finite difference
        Returns :
                matrix : np matrix
        """
        N = self.npts
        D2 = np.zeros(shape=(N,N))
        for i in range(1,N-1):
            D2[i,i]   = -2 
            D2[i,i-1] = 1 
            D2[i,i+1] = 1 
        
        D2[0,0]  = D2[N-1,N-1] = -2
        D2[0,1]  = D2[N-1,N-2] = 1
        
        D2 /= self.delta**2 
        return D2

    
    def centered_1st_derivative(self):
        """
        Desc : Generates a tri-diagonal matrix for 1st order 
               derivative using a finite difference
        Returns :
                matrix : np matrix
        """
        N = self.npts + 2
        D1 = np.zeros(shape=(N,N))
        for i in range(1,N-1):

            D1[i,i]   = 0 
            D1[i,i-1] = -1 
            D1[i,i+1] = 1 
        
        #One side boundary derivative
        D1[0,0]  = D1[N-1,N-1] = 0
        D1[0,1]     = 4
        D1[0,2]     = -1 
        D1[N-1,N-2] = -4
        D1[N-1,N-3] = 1
        
        D1 /= 2*self.delta 
        return D1
    
    
    def centered_1st_derivative_fs(self):
        """
        Desc : Generates a tri-diagonal matrix for 1st order 
               derivative using a finite difference
        Returns :
                matrix : np matrix
        """
        N = self.npts + 2
        D1 = np.zeros(shape=(N,N))
        for i in range(1,N-1):
            D1[i,i]   = 0 
            D1[i,i-1] = -1 
            D1[i,i+1] = 1 
        
        #One side boundary derivative
        D1[N-1,N-1] = 0
        D1[N-1,N-2] = -4
        D1[N-1,N-3] = 1
        
        D1 /= 2*self.delta 
        D1[0,0]  = -1/self.delta
        return D1
        
    
    def evaluate_mode(self,mode_number,z):
        if self.mode_funcs:
            return self.mode_funcs[mode_number](z)
            
        else:
            self.interpolate_modes()
            return  self.mode_funcs[mode_number](z)
    
    
    def interpolate_modes(self):
        self.mode_funcs = []
        for n in range(self.modes.shape[1]):
            import scipy
            fun = scipy.interpolate.interp1d(self.domain,self.modes[:,n])
            self.mode_funcs.append(fun)
    
    
class InternalWaveModes(ModeSolver):
    def __init__(self,N2,Z,frequency,
                           latitude=30,
                           scalar_gradient=None,
                           free_surface=False,
                           puvmodes=False):
        
        self.fo = coriolis(latitude)
        self.frequency = frequency
        if ( self.frequency < self.fo ):
            print(frequency,self.fo)
            print("Frequency needs to be greater than coriolios")
            raise
        
        elif( self.frequency > 2*np.sqrt(max(N2)) ):
            print(frequency,np.sqrt(max(N2)))
            print("Frequency cannot exceed max stratfication")
            raise
        
        elif( self.frequency > np.sqrt(max(N2)) ):
            print("Halving frequency for dispersion's sake")
            print(self.frequency," => ", self.frequency/2)
            frequency /= 2
            
        F = frequency**2 * np.ones(len(Z))    
        LH = np.diag(N2-F)
        super().__init__(Z,LH,boundary=[0,0],free_surface=free_surface)
        
        self.hwavenumbers, self.modes = self.solve(frequency)
        
        if puvmodes:
            if free_surface:
                D1 = self.centered_1st_derivative_fs()
            else:
                D1 = self.centered_1st_derivative()
                
            self.modes = D1 @ self.modes 
        
        if scalar_gradient is not None:
            #Multiple rows by each element in gradient
            self.modes = (self.modes.T * scalar_gradient).T
    
    def solve(self,frequency=None):
        eigval, eigvec = super().solve()
        frequency = self.frequency if frequency is None else frequency
        eigval = np.sqrt((frequency**2 - self.fo**2)/
                         ( (2*np.pi)**2 *eigval))
        return eigval,eigvec
    
    
    def hwavenumber(self,mode_number):
        return self.hwavenumbers[mode_number].real
    

class QuasiGeostrophicModes(ModeSolver):
    def __init__(self,N2,Z,frequency=0,scalar_gradient=None,free_surface=False,latitude=30):
        LH = np.diag(N2)
        self.fo = coriolis(latitude)
        self.frequency=frequency
        super().__init__(Z,LH,free_surface=free_surface)

        self.hwavenumbers, self.modes = self.solve()
        
        if scalar_gradient is not None:
            #Multiple rows by each element in gradient
            self.modes = (self.modes.T * scalar_gradient).T
    
    
    def solve(self):
        eigval, eigvec = super().solve()
        eigval = self.fo/np.sqrt(eigval)
        return eigval,eigvec

    
    def hwavenumber(self,mode_number):
        return self.hwavenumbers[mode_number]
    
    
def coriolis(latitude):
    #rotation rate of earth in cps
    omega = 7.292115*1e-5/(2*np.pi) 
    f = 2*omega*np.sin(np.pi*latitude/180)
    return f