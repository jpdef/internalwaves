import numpy as np 

def tappered_least_square(H,RINV=None,QINV=None):    
    RINV = RINV if RINV is not None else np.diag(np.ones(H.shape[0]))
    QINV = QINV if QINV is not None else np.diag(np.ones(H.shape[1]))
    
    HDAG  = H.T @ RINV @ H + QINV
    HDAG  = np.linalg.inv(HDAG) @ H.T @ RINV
    
    return HDAG