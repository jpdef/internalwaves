import numpy as np 

def tappered_least_square(H,RINV=None,QINV=None):
    QINV = QINV if QINV is not None else np.diag(np.ones(H.shape[1]))
    
    if RINV is not None:
        M = np.sqrt(RINV) * H
        HDAG = np.linalg.inv(M.T @ M + np.diag(QINV) ) @ M.T * np.sqrt(RINV)
        
    else:
        HDAG = np.linalg.inv(H.T @ H + np.diag(QINV) ) @ H.T
        
        
    RINV = RINV if RINV is not None else np.diag(np.ones(H.shape[0]))
    
    return HDAG

