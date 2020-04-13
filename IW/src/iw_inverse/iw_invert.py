#Desc : Library for evaluating parameters of observed internal wave field
#Auth : J. DeFilippis
#Date : 10-28-2019

import pandas as pd
import numpy as np
import itertools

import sys
sys.path.append('../src')

import pseudo_inverse as pinv

from iw_modes import InternalWaveModes


class InternalWaveInversion:
    """
    Desc: Class to invert for parameters frequency,mode, and wavelength
    from oceanographic data using a linear internal wave model
    """
    
    def __init__(self,ds,freqs,modes,thetas,depth,N2=np.array([]),tc='dz',mean=True):
        """
        Desc: Sets inputs to field values. Creates a dataframe of parameters
        to invert for. Creates a model matrix from data coordinates x,y,z,t.
        """
        self.set_attributes(ds,freqs,modes,thetas,depth,N2,tc)
        self.ps = self.generate_param_space()
        self.H  = self.model_matrix()
        
        #Fit a mean if the frequency is zero
        if mean:
            new_row = pd.DataFrame({"freqs" : [0], "modes" : [0], "hwvn" : [0], "theta" : [0]})
            self.ps = pd.concat([new_row, self.ps[:]]).reset_index(drop = True) 
            
            mid = int( np.floor((self.H.shape[1])/2) ) + 1
            self.H = np.insert(self.H,0,np.ones(self.H.shape[0]).T,axis=1)
            self.H = np.insert(self.H,mid,np.ones(self.H.shape[0]).T,axis=1)
    
    
    def set_attributes(self,ds,freqs,modes,thetas,depth,N2,tc):
        self.freqs  = freqs
        self.modes  = modes
        self.N2     = N2
        self.depth  = depth
        self.ds     = ds
        self.tc     = tc
        self.thetas = thetas
    
    
    def generate_param_space(self):
        #Create combination of frequency,modes, and wave direction
        parr = np.array(list(itertools.product(*[self.freqs,self.modes,self.thetas])),
                         dtype=[('f',float),('m',int),('theta',float)])
        
        #Generate wavenumbers from disperion relation
        self.iwmodes = [InternalWaveModes(self.depth,self.N2,freq=f) for f in self.freqs]
    
        #Map wavenumbers to mode,frequency
        hwvn = [self.iwmodes[ self.freqs.index(p['f']) ].get_hwavenumber(p['m']) for p in parr] 
        
                        
        return pd.DataFrame({'freqs' : parr['f'], 
                             'modes' : parr['m'], 
                             'hwvn'  : hwvn,
                             'theta' : parr['theta']
                            })
        

    def model_matrix(self):
        vmat = self.vertical_matrix()
        hmat = self.horizontal_matrix()
        
        return np.multiply(vmat,hmat)
        
        
    def horizontal_matrix(self):
        """
        Desc: Create horizontal portion of matrix with sin and cos basis functions. 
        ie :// a sin(k . x  - ft) + b cos( k . x - ft)
        """
        #Make arguement
        rad =  self.ps['theta'] * (np.pi/180)
        kx = 2*np.pi*np.outer(self.ds['x'],self.ps['hwvn']*np.cos(rad)) 
        ky = 2*np.pi*np.outer(self.ds['y'],self.ps['hwvn']*np.sin(rad)) 
        ft = 2*np.pi*np.outer(self.ds['time'], self.ps['freqs'])
        
        #Compute for sin and cos parts
        psi  = np.hstack([np.sin(kx + ky - ft),np.cos(kx + ky - ft)])
        psi = psi.astype('float')
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
        phi = phi.T
        phi = phi.astype(float)
        return phi


    def pinvert(self):
        if 'RINV' in self.ds:
            RINV = np.diag(self.ds['RINV'])  
        else:
            RINV = np.diag(np.ones(self.H.shape[0]))
        
        if 'QINV' in self.ps:
            q = np.concatenate([self.ps['QINV'],self.ps['QINV']])
            QINV = np.diag(q)  
        else:
            QINV = np.diag(np.ones(self.H.shape[1]))
        
        self.hdag = pinv.tappered_least_square(self.H,RINV=RINV,QINV=QINV)
            
        amps = self.hdag @ self.ds[self.tc] 
        div  = int(np.floor(len(amps)/2))
        self.ps['a'] = amps[0:div]
        self.ps['b'] = amps[div:]



        
        
        
