from scipy import interpolate
import numpy as np 
from tqdm import tqdm


class ScalarField:
    def __init__(self,scalar_prof,depth):
        self.init_prof = scalar_prof;
        self.depth = depth
        self.depth_map = interpolate.interp1d(depth,scalar_prof)
        
        
    def advect(self,displacements):
        """
        Displacements are a 2d array axis : 0 =  time , 1=depth
        Returns 2d slice of scalar field
        """
        #Iterate over time
        scalar_prof = self.init_prof
        scalar_field = np.zeros(displacements.shape)
        for i in tqdm(range(displacements.shape[0])):
            #Update Z
            delta_z = displacements[i,:]
            zprime = self.depth - delta_z
            if self.unstable(zprime):
                print('Inversion Occurred')
                raise
            
            #Update DMAP
            dmap  = interpolate.interp1d(zprime,scalar_prof)
            
            #Map back to regular grid
            scalar_prof = dmap(self.depth)
            scalar_field[i,:]  = scalar_prof
        
        return scalar_field
    
    
    def advect_bruce(self,displacements):
        """
        Displacements are a 2d array axis : 0 =  time , 1=depth
        Returns 2d slice of scalar field
        """
        #Iterate over time
        scalar_field = np.zeros(displacements.shape)
        for i in tqdm(range(displacements.shape[0])):
            #Update Z
            delta_z = displacements[i,:]
            zprime = self.depth - delta_z
            if self.unstable(zprime):
                print('Inversion Occurred')
                raise
            
            #Update DMAP
            dmap  = interpolate.interp1d(zprime, self.init_prof)
            
            #Map back to regular grid
            scalar_field[i,:] = dmap(self.depth)
        
        return scalar_field
        
        
    def unstable(self,z):
        zprev = -1
        for zv in z:
            if zv < zprev:
                print("Unstable : ",zv,zprev)
                return True
            zprev = zv
            
        return False