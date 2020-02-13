#Desc : Library for evaluating parameters of observed internal wave field
#Auth : J. DeFilippis
#Date : 10-28-2019

import pandas as pd
import numpy as np

import sys
sys.path.append('../src')

import pseudo_inverse as pinv

from iw_modes import InternalWaveModes


class InternalWaveInversion:
    """
    Desc: Class to invert for parameters frequency,mode, and wavelength
    from oceanographic data using a linear internal wave model
    """
    
    def __init__(self,ds,freqs,modes,depth,N2=np.array([]),tc='dz'):
        """
        Desc: Sets inputs to field values. Greats a dataframe of parameters
        to invert for. Creates a model matrix from data coordinates x,y,z,t.
        """
        self.set_attributes(ds,freqs,modes,depth,N2,tc)
        self.ps = self.generate_param_space()
        self.H  = self.model_matrix()
     
    
    def set_attributes(self,ds,freqs,modes,depth,N2,tc):
        self.freqs  = freqs
        self.modes  = modes
        self.N2     = N2
        self.depth  = depth
        self.ds     = ds
        self.tc     = tc
    
    
    def generate_param_space(self):
        #Make columns a of frequencies and modes in correct index order
        freqs = np.repeat(self.freqs,len(self.modes))
        modes = np.tile(self.modes,len(self.freqs))
        
        #Generate column of wavenumbers
        self.iwmodes= [InternalWaveModes(self.depth,self.N2,freq=f) for f in self.freqs]
        hwvn = np.array([ self.iwmodes[i].get_hwavenumber(self.modes) 
                         for i in range(len(self.freqs))]).flatten()
        
        return pd.DataFrame({'freqs' : freqs, 'modes' : modes , 'hwvn' : hwvn})
        

    def model_matrix(self):
        vmat = self.vertical_matrix()
        hmat = self.horizontal_matrix()
        
        return np.multiply(vmat,hmat)
        
        
    def horizontal_matrix(self):
        """
        Desc: Create horizontal portion of matrix with sin and cos basis functions. 
        ie :// a sin(k . x - ft) + b cos( k . x - ft)
        """
        #Make arguement
        kx = 2*np.pi*np.outer(self.ds['x'],self.ps['hwvn']) 
        ft = 2*np.pi*np.outer(self.ds['time'], self.ps['freqs'])
        
        #Compute for sin and cos parts
        psi  = np.hstack([np.sin(kx + ft),np.cos(kx+ft)])
        self.psi = psi
        return psi
    
    
    def vertical_matrix(self):
        """
        Desc: Create  portion of matrix with vertical basis functions. 
        ie :// phi_nm (z)
        """
        cols = []
        
        #Get index for frequency values
        def find_nearest(array,value):
            idx = (np.abs(array-value)).argmin()
            return idx
        
        #Iterate over parameter rows to get correct modal function
        #put them in a list
        for index,row in self.ps.iterrows():
            f = row['freqs'] 
            m = int(row['modes'].real)
            i = find_nearest(self.freqs,f)
            col = self.iwmodes[i].evaluate_mode(m,self.ds['z'])
            cols.append(col)
            
            
        """"
        Make a matrix rows of phi(z) then duplicate for cosine part
        and transpose so parameters are in column space and data
         is in row space 
        """
        phi = np.vstack(cols)
        phi = np.vstack([phi,phi])
        self.phi = phi.T
        return phi.T


    def pinvert(self):
        self.hdag = pinv.tappered_least_square(self.H)
        amps = self.hdag @ self.ds[self.tc] 
        div  = int(np.floor(len(amps)/2))
        self.ps['a'] = amps[0:div]
        self.ps['b'] = amps[div:]


def read_data_dir(path):
    for f in get_file_list(path):
        pass

def get_file_list(path,fpat='run*'):
    file + glob.glob(os.path.join(path,fpat))
    return sorted (files)

