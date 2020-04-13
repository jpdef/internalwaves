import numpy as np

#Constants
b  = 0 
No = 0
N  = 0
E  = 6.3e-5
f  = 7.3e-5
jstar = 3

def B(freq):
    return (2*f/(np.pi*freq**3))*np.sqrt(freq**2-f**2) 

def H(j):
    M = 0.5*jstar**2/(np.pi*jstar-1)
    return M/(j**2 + jstar**2) 

def gm_spectral_dist(freq,j):
    return B(freq)*H(j)

def generate_gm_amp(freq,j):
    var = gm_spectral_dist(freq,j)
    a = np.random.normal(0,np.sqrt(var),1)
    return(a)

def generate_gm_pm(freq,j):
    a = generate_gm_amp(freq,j)
    b = generate_gm_amp(freq,j)
    phase = np.arctan(b/a)
    mag   = np.sqrt(b**2 + a**2)
    return(mag,phase)


#Munk sound depth profile
def sound_prof_munk( z ):
    cbar = 1450
    epi = 0.00737
    return cbar*(1 + epi*(zbar(z)- 1 + np.exp(-zbar(z))))

def sound_grad_munk(z):
    cbar = 1500
    epi = 0.00737
    return cbar*(2/1300)*epi*( 1 - np.exp(-zbar(z)) )

#Average depth described by munk profile
def zbar (z):
    return 2*(z-1300)/1300

@np.vectorize
def haversine(lon1, lat1, lon2, lat2):
    from math import radians, cos, sin, asin, sqrt
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    def torad(angle):
        return 180*angle/np.pi
    
    lon = torad(lon1) 
    lat = torad(lat1)
    lon = torad(lon2) 
    lat = torad(lat2)

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r