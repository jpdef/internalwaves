#SCALAR_FIELD
#Desc : Two dimensional representation of a scalar field
#Auth : J. DeFilippis
#Date : 7-29-2019

import numpy as np

class ScalarField:
    def __init__(self,profile,npoints):
        if profile.size != npoints[0] :
            print("Dimensions of depth profile and field points not the same")
        
        self.size = npoints[0]

        #Create field matrix (copies depth profile n times)
        self.field = np.array([profile,]*npoints[1]).transpose()
        self.base = np.array([profile,]*npoints[1]).transpose()

        #Create gradient matrix
        self.gradient_field = self.make_gradient() 

    def distort(self,displacement_field):
        #Update Field
        self.field = self.base + np.multiply(displacement_field,self.gradient_field)
        
    def make_gradient(self):
        #Compute gradient
        grad =  np.gradient(self.field)
        
        # Return column wise gradient
        return grad[0] 
