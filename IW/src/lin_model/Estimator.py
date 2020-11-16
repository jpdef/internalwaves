#Desc : Library for linear estimation
#Auth : J. DeFilippis
#Date : 7-25-2020

import numpy as np 
import LinearModel 
from tqdm import tqdm
from functools import partial

class Estimator:
    """
    Desc:
    A class to estimate parameters of a linear model
    @ linearmodel (LinearModel) model to to estimate parameters
    @ Q (array) prior uncertianity of model parameters (diagonal of covariance)
    """
    def __init__(self,linearmodel,Q=None):
        self.linearmodel = linearmodel
        
        if Q is None:
            self.Q = np.ones(linearmodel.columns)
        elif (any(Q) > 0):
            self.Q = Q
        else:
            print("Q must be Q > 0 ")
            raise
        
        
    def estimate(self,coordinates,var='dz'):
        """
        Estimates the parameters of linear model by multiplying
        adjoint onto observations.
        @param coordinates
               list of tuple with x,y,z,t and observation variable
        @param var
               name of the observation to fit in tuple list
        """
        H = self.linearmodel.model_matrix(coordinates)
        Ht = self.adjoint(H,None)
        return Ht@coordinates[var]
   
            
    def adjoint(self,H,R):
        """
        Computes the adjoint of model H with weighted
        matrices R and Q.
        Ht = (H.T R^-1 H + Q^ -1)^-1 H.T
        M = (sqrt(R) * H) (for efficiency)
        """
        if R is None:
            QI = np.diag(1/self.Q)
            Ht = np.linalg.inv(H.T @ H + QI ) @ H.T
            return Ht
            
        r = np.diag(1/R)
        M = ( np.sqrt(r)*H )
        Ht = np.linalg.inv(M.T @ M + QI ) @ M.T * r
            
        return Ht

    
    def recursive_estimate(self,coordinate_list,var='dz'):
        """
        Computes a recursive least squares (RLS) esimate on data. This is
        suggested for large (N>100) datasets. Reduces size of matrix inversion.
        @param coordinates_list
               list of list of tuple with x,y,z,t and observation variable
               first list segements the chunks of data to invert. 
        @param var
               name of the observation to fit in tuple list
        
        """
        parameter_estimate = np.random.rand(self.linearmodel.columns)
        Q = np.diag(self.Q)
        for i,coordinates in tqdm(enumerate(coordinate_list)):
            H = self.linearmodel.model_matrix(coordinates)
            estimate = H @ parameter_estimate
            innovation = coordinates[var] - estimate
            
            if 'R' in coordinates.dtype.names:
                R = np.diag(coordinates['R'].flatten())
            else: 
                R = np.identity(len(coordinates))
                
            #Correct & Update
            INV = np.linalg.inv( H @ Q @ H.T + R )
            K = Q @ H.T @ INV
            Q = (Q - K @ H @ Q)
            parameter_estimate = parameter_estimate + K @ innovation
        
        return [parameter_estimate,Q]
    
    
    def recursive_posterior_uncertianity(self,coordinate_list):
        """
        Computes a posterior uncertianity of least squares estimate. This is
        suggested for large (N>100) datasets. Reduces size of matrix inversion.
        @param coordinates_list
               list of list of tuple with x,y,z,t and observation variable
               first list segements the chunks of data to invert.         
        """
        Q = np.diag(self.Q)
        for i,coordinates in tqdm(enumerate(coordinate_list)):
            H = self.linearmodel.model_matrix(coordinates)
            
            if 'R' in coordinates.dtype.names:
                R = np.diag(coordinates['R'].flatten())
            else: 
                R = np.identity(len(coordinates))
                
            #Correct & Update
            INV = np.linalg.inv( H @ Q @ H.T + R )
            K = Q @ H.T @ INV
            Q = (Q - K @ H @ Q)
            #print(np.trace(Q))
        
        return Q
   

    def posterior_uncertianity(self,coordinates):
        """
        Computes a posterior uncertianity of least squares estimate. 
        @param coordinates
               list of tuple with x,y,z,t and observation variable
        """
        Q = np.diag(self.Q)
        H = self.linearmodel.model_matrix(coordinates.flatten())
            
        if 'R' in coordinates.dtype.names:
            R = np.diag(coordinates['R'].flatten())
        else: 
            R = np.identity(len(coordinates))
                
        INV = np.linalg.inv( H @ Q @ H.T + R )
        K = Q @ H.T @ INV
        Q = (Q - K @ H @ Q)
        
        return Q

    
    def error_map(self,coordinates,Q=None):
        """
        Computes an error map. Model space uncertianity transformed with
        the model into data space
        @param coordinates
               list of tuple with x,y,z,t and observation variable
        """
        H = self.linearmodel.model_matrix(coordinates)
        if Q is None:
            return  np.copy(np.diagonal(H @ np.diag(self.Q) @ H.T))
        else:
            return  np.copy( np.diagonal(H @ Q @ H.T)  )