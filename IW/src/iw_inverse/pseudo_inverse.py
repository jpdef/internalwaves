import numpy as np 

def tappered_least_square(H) :
    RINV = np.diag(np.ones(H.shape[0]))
    QINV = np.diag(np.ones(H.shape[1]))
    
    HDAG  = H.T @ RINV @ H + QINV
    HDAG  = np.linalg.inv(HDAG) @ H.T @ RINV
    
    return HDAG