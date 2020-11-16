#Desc : Library for Internal Wave Spectrum Modeling
#Auth : J. DeFilippis
#Date : 9-28-2020
import sys
import numpy as np 
import scipy
sys.path.append('../src/misc')
sys.path.append('../src/lin_model')
from iw_misc import coriolis
from lin_model import LinearModel 
from lin_model import ModeSolver
from tqdm import tqdm
from functools import partial
tocph = 3600/(2*np.pi)

"""
Garrett Munk Meta Parameters 
"""
class GM81Params:
    """
    Meta parameters for Garrett-Munk model 
        E  energy scale constant
        b  depth of thermocline
        No maximum (surface) stratification
    """
    def __init__(self):
        self.E = 6.3e-5
        self.b = 1.3e3
        self.No= 3/(tocph)

        
class CustomParams:
    def __init__(self,E,N2,Z):
        self.E = E
        self.No, self.b = stratification_fit(N2,Z)
    
    def stratification_fit(self,N2,Z):
        def func(x, a, b, c):
            return a * np.exp(-b*x) + c
        N = np.sqrt(N2)
        popt,pcov = curve_fit(func,Z,N,p0=[np.max(N),Z[len(Z)//2],0])
        return popt[0],popt[1]

    
    
    
class ModalSpectralModel:
    """
    Internal Wave Power Spectrum Using EigenModes.
        F_d = B(omega)*H(j)*R(omega,j)
        @param N2       : stratification  in [rad/s] nx1 
        @param Z        : depth in [m] nx1
        @param latitude : local line of latitude in [deg]
        @param jstar    : parameters for IW mode power distribution
        @param nmodes   : number of modes to include in spectrum
    """
    def __init__(self,N2,Z,latitude=30,jstar=3,nmodes=10):
        #Assign fields
        self.N2 = N2
        self.Z  = Z 
        self.jstar  = jstar
        self.nmodes = nmodes
        self.C  = 1/np.sum([self.H(j) for j in range(1,nmodes+1)])
        self.latitude=latitude
        self.f = (2*np.pi)*coriolis(latitude)
        self.A = 1
        self.energy_density  = 3800 # J/m^2
        self.mean_sw_density = 1025
                
    
    def energy_normal(self,iwm):
        """
        Normalizes displacement modes to have energy that corresponds with
        Garret-Munk
        """
        a = np.zeros(self.nmodes)
        for i in range(0,self.nmodes):
            integral = scipy.integrate.trapz(self.N2*(iwm.modes[:,i])**2,self.Z)
            a[i] = (2*self.energy_density/(self.mean_sw_density*integral) )
            #a[i] = 1/integral
        return np.sqrt(a)
    
    
    def evaluate_modal_amplitudes(self,omega):
        """
        Computes modes and modal amplitudes per frequency
        """
        iwm = ModeSolver.InternalWaveModes(self.N2,self.Z,omega,latitude=self.latitude)
        amplitudes = self.energy_normal(iwm)
        return iwm,amplitudes
    
    
    def H(self,j):
        """
        Garrett-Munk Mode Distribution
        """
        return 1/(j**2 + self.jstar**2)
    
    
    def R(self,omega,j):
        """
        Rob Pinkel correction to Garret Munk 
        """
        threshold = 1/tocph # 1cph
        sig = 1 / (1 + np.exp(-100*(omega-threshold)*tocph))
        rolloff = (1-sig) + sig*np.exp(-(j-1)*(omega-threshold)*tocph)
        return rolloff
    
    
    def Comega(self,omega):
        return 1/np.sum([self.H(j)*self.R(omega,j) for j in range(1,self.nmodes+1)])
    
    #Frequency distribution
    def B(self,omega):
        """
        Garrett-Munk Frequency Distribution
        """
        return 2*self.f/(np.pi*omega*np.sqrt(omega**2-self.f**2)) 
    
    
    #Internal Wave Displacement Power Spectral Density
    def FD(self,omega,depth,mode=None,): 
        power = np.zeros(len(omega))
        #Case 1: Mode is specified
        if mode is not None:
            for oi,o in enumerate(omega):
                #Compute mode and mode amplitude evaluate those at a specific depth
                vertical_modes, modal_amplitudes = self.evaluate_modal_amplitudes(o)
                AofZ = (modal_amplitudes[mode]*vertical_modes.evaluate_mode(mode,depth))**2
                
                #Compute GM power density and multiply by modal amplitude
                P =  self.B(o)* self.H(mode)*self.C
                power[oi] =  AofZ*( (o**2 - self.f**2)/o**2 ) * P
        
        #Case 2: Mode is not specified
        else:
            for oi,o in enumerate(omega):
                #Compute mode and mode amplitude evaluate those at a specific depth
                vertical_modes, modal_amplitudes = self.evaluate_modal_amplitudes(o)
                AofZ = (modal_amplitudes[mode]*vertical_modes.evaluate_mode(mode,depth))**2
                
                #Compute GM power density and multiply by the sum of modal amplitudes
                P =  self.B(o)*self.C
                SumAofZ = np.sum([self.H(m+1)* (modal_amplitudes[m]*vertical_modes.evaluate_mode(m,depth))**2 for m in range(0,self.nmodes)])
                power[oi] =  SumAofZ*( (o**2 - self.f**2)/o**2 ) * P
            
        return  self.A*power
    
    
    #Internal Wave Displacement Power Spectral Density with Mode Roll Off 
    def FDRP(self,omega,depth=None,mode=None,): 
        if depth is not None:
            power = np.zeros(len(omega))
            for oi,o in enumerate(omega):
                iwm = ModeSolver.InternalWaveModes(self.N2,self.Z,o,latitude=self.latitude)
                amplitudes = self.energy_normal(iwm)
                Co  = self.Comega(o)*self.C
                P =  self.B(o)* self.H(mode)*self.R(o,mode)*Co
                asq = (amplitudes[mode]*iwm.evaluate_mode(mode,depth))**2
                power[oi] =  asq*( (o**2 - self.f**2)/o**2 ) * P
        else:
            for oi,o in enumerate(omega):
                iwm = ModeSolver.InternalWaveModes(self.N2,self.Z,o,latitude=self.latitude)
                amplitudes = self.energy_normal(iwm)
                Co  = self.Comega(o)
                P =  self.B(omega)* self.H(mode)*self.R(o,mode)*Co
                asq = np.sum([self.H(m+1)*self.R(o,m) (amplitudes[m]*iwm.evaluate_mode(m,depth))**2 for m in range(0,self.nmodes)])
                power[oi] =  asq*( (o**2 - self.f**2)/o**2 ) * P
            
        return  self.A*power
 
    def FDOM(self,omega,mode):
        P = self.B(omega)*self.H(mode)*self.C
        return P*( (omega**2 - self.f**2)/omega**2 ) 


    def FK(self,depth,omegas,modes,angles,sg=None):
        alphas = np.zeros(shape=(len(omegas),len(modes),len(angles)))
        power  = np.zeros(shape=(len(omegas),len(modes),len(angles)))
        zi = np.argmin( abs( self.Z - depth ))
        #Iterate over frequency
        for oi,o in enumerate(omegas):
            iwm = ModeSolver.InternalWaveModes(self.N2,self.Z,o) 
            amplitudes = self.energy_normal(iwm)
            Co = self.Comega(o)
            
            #Iterate over mode
            for mi,m in enumerate(modes):
                kmag = iwm.hwavenumber(m) 
                S = 1 if sg is None else sg[zi] 
                mode_amp = (S*amplitudes[m]*iwm.modes[zi,m])**2
                P = mode_amp*self.B(o)*self.H(m+1)*self.R(o,m+1)*Co
                
                #Iterate over angle
                for ai,a in enumerate(angles):
                    alphas[oi,mi,ai] = kmag*np.cos(a*np.pi/180)
                    power[oi,mi,ai]  = P/len(angles)
        
        return alphas.flatten(),power.flatten()
    
    
    def wavenumber_space(self,iwm):
        kx = np.zeros(shape=(len(freqs),len(angles),len(modes)))
        ky = np.zeros(shape=(len(freqs),len(angles),len(modes)))
        for j,a in enumerate(angles):
            kmag = np.array([iwm.hwavenumber(m) for m in range(self.nmodes)])
            kx[i,j,:] = np.cos(a*np.pi/180)*kmag
        return kx
    
    
    def wavenumber_power(self,iwm):
        power = np.zeros(shape=( len(freqs),len(angles),len(modes)) )
        amplitudes = self.energy_normal(iwm)
        for m,mode in enumerate(modes):
            #Power at depth zi
            power_per_depth = (amplitudes[m]*iwm.modes[zi,m])**2
            psd = B(omega)*R(omega,mode)*Co
            for j,a in enumerate(angles):
                power[i,j,m] = (1/tocpkm)*power_per_depth*psd/len(angles)
        return power
    
        
    
class WKBJSpectralModel:
    """
    Internal Wave Power Spectrum Using WKBJ.
        F_d = B(omega)*H(j)*R(omega,j)
        @param N2       : stratification  in [rad/s] nx1 
        @param Z        : depth in [m] nx1
        @param latitude : local line of latitude in [deg]
        @param jstar    : parameters for IW mode power distribution
        @param nmodes   : number of modes to include in spectrum
    """
    def __init__(self,smparams,N,latitude,jstar,nmodes=10): 
        self.E  = smparams.E
        self.b  = smparams.b
        self.No = smparams.No
        self.N  = N
        self.nmodes = nmodes
        self.jstar = jstar
        self.f  = (2*np.pi)*coriolis(latitude) 
        self.C  = 1/np.sum([self.H(j) for j in range(1,nmodes+1)])
        self.A = self.E*(self.b**2)/(2*self.No/self.N)
        print(self.f,self.C,self.A)
    
    #Modal distribution
    def H(self,j):
        return 1/(j**2 + self.jstar**2)
    
    #Frequency distribution
    def B(self,omega):
        return 2*self.f/(np.pi*omega*np.sqrt(omega**2-self.f**2)) 
    
    
    #Internal Wave Displacement Power
    def FD(self,omega,alpha=None,mode=None,):
        if alpha is not None:
            ktom = (self.b/np.pi)*np.sqrt( (self.No**2-omega**2)/(omega**2-self.f**2) )
            mode = np.pi*ktom * k 
            P = (1/ktom)*self.B(omega)*self.H(mode)*self.C
            
        elif mode is not None:
            P =  self.B(omega)* self.H(mode)*self.C
        
        else:
            P =  self.B(omega) 
            
        return  self.A * ( (omega**2 - self.f**2)/omega**2 )*P
                 
        
class GM81(WKBJSpectralModel):
    def __init__(self,N,smparams=GM81Params(),latitude=30,jstar=2,nmodes=10):
        super().__init__(smparams,N,latitude,jstar,nmodes=nmodes)
     
    
class RP20(WKBJSpectralModel):
    """
    Rob pinkel internal wave spectral model
    """
    def __init__(self,N,smparams=GM81Params(),latitude=30,jstar=2,nmodes=10):
        super().__init__(smparams,N,latitude,jstar,nmodes=nmodes)
    
    def R(self,omega,j):
        threshold = 1/tocph # 1cph
        sig = 1 / (1 + np.exp(-100*(omega-threshold)*tocph))
        rolloff = (1-sig) + sig*np.exp(-(j-1)*(omega-threshold)*tocph)
        return rolloff
    
    def Comega(self,omega):
        return 1/np.sum([self.H(j)*self.R(omega,j) for j in range(1,self.nmodes+1)])
   

    #Internal Wave Displacement Power
    def FD(self,omega,alpha=None,mode=None,):
        Co  = [(self.Comega(o)) for o in omega]
        if alpha is not None:
            ktom = (self.b/np.pi)*np.sqrt( (self.No**2-omega**2)/(omega**2-self.f**2) )
            mode = np.pi*ktom * k 
            P = (1/ktom)*self.B(omega)*self.H(mode)*self.C
            
        elif mode is not None:
            P =  self.B(omega)* self.H(mode)*self.R(omega,mode)*Co
        
        else:
            P =  self.B(omega) 
            
        return self.A * ( (omega**2 - self.f**2)/omega**2 ) *P
    
    
        
        
