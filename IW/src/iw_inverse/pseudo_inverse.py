import numpy as np 

def tappered_least_square(H,RINV=None,QINV=None):
    QINV = np.diag(QINV) if QINV is not None else np.diag(np.ones(H.shape[1]))

    if RINV is not None:
        M = (np.array(H.T) *  np.sqrt(RINV)).T
        HDAG = np.linalg.inv(M.T @ M + QINV ) @ M.T * np.sqrt(RINV)
        
    else:
        HDAG = np.linalg.inv(H.T @ H + QINV ) @ H.T
        
        
    return HDAG


def left_inverse(H,RINV=None,QINV=None):
    QINV = np.diag(QINV) if QINV is not None else np.diag(np.ones(H.shape[1]))

    if RINV is not None:
        M = (np.array(H.T) *  np.sqrt(RINV)).T
        HDAG = np.linalg.inv(M.T @ M + QINV ) @ M.T * np.sqrt(RINV)
        
    else:
        HDAG = np.linalg.inv(H.T @ H + QINV ) @ H.T
        
        
    return HDAG



def right_inverse(H,RINV=None,QINV=None):
    RINV = np.diag(RINV) if RINV is not None else np.diag(np.ones(H.shape[0]))

    if True:
        M = (np.array(H) @  np.diag(np.sqrt(QINV)) ).T
        INV= np.linalg.inv(M.T @ M + RINV ) 
        HDAG =   (np.diag(np.sqrt(QINV) )@ M) @ INV
        
    else:
        HDAG = H.T @ np.linalg.inv(H @ H.T + RINV )
        
        
    return HDAG