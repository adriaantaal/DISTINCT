import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import scipy as sp
from scipy import io as spio
from scipy import optimize as spop

## The code for distinct as debugged and simulated in the 
def DSv7(yf2, D, Ek, lambdaval, H, g, intSolver, maxIter, maxIterInt, addSize):
    print("Running DISTINCT ....")
    M,N = D.shape
    K,N = Ek.shape
    
    #initialize
    incrowd = np.zeros(1, dtype=np.int64)
    oldcrowd = np.zeros(1, dtype=np.int64)
    J = np.zeros((maxIter,1), dtype=np.int64)
    xhatFull = np.zeros((N, 1), dtype=float)
    u = g
    
    #start for loop
    for ii in range(0,maxIter-1):
        us = np.flipud(np.sort(u))
        us[addSize+1:N] = 0
        indus = np.flipud(np.argsort(u))
        toAdd = indus[us>0]
        
        #detect new crowd
        if ii == 0:
            incrowd = toAdd
        else:
            ictoAdd = np.concatenate([incrowd,toAdd])
            _, idx = np.unique(ictoAdd, return_index=True)
            incrowd = ictoAdd[np.sort(idx)]
            
            
        #run the solver
        if intSolver == 0:
            x0 = np.zeros(incrowd.shape, dtype=float)
            
            if ii > 0:
                oldcrowdlen = oldcrowd.shape
                x0[1:oldcrowdlen[0]] = xhat
            
            Hsmall = H[incrowd,:]
            Hsmall= Hsmall[:,incrowd]
            xhat = spop.nnls(Hsmall, g[incrowd], maxiter=maxIterInt)
            xhat = xhat[0]
            # np.asarray(xhat,dtype=float)
            
        # remove zero entries        
        
        # update u
        u = g - Hsmall@xhat
            
        #loop back
        oldcrowd = incrowd
        
        J[ii] = yf2 + np.transpose(xhat)@(Hsmall@xhat) + lambdaval*sum(abs(xhat))
        
        #find if J converged
        if ii>1:
            if (abs(J[ii]-J[ii-1])/J[ii-1] < 1e-3) || (abs(J[ii]/yf2) < 1e-4):
                break
        

    xhatFull[incrowd] = xhat
        
    return xhatFull, J
    
def DSv7TB():
    
    # M = int(1024)
    # K = int(1024)
    # N = int(6468)
    # D = np.random.randn(M,N)
    # Ek = np.random.randn(K,N)
    
    #get matrices from .mat file
    mat_contents = spio.loadmat("C:\DSmats")
    print(mat_contents.keys())
    D = mat_contents.get("D")
    Ek = mat_contents.get("Ek")
    X = mat_contents.get("X")
    Y = mat_contents.get("Y")
    Z = mat_contents.get("Z")
    
    #derive from matrices
    M,N = D.shape
    K,N = Ek.shape
    DTD = np.transpose(D)@D
    ETE = np.transpose(Ek)@Ek
    H = DTD*ETE #hadamard product
    
    #generate data
    lambdadiv = float(10)

    y = np.random.randn(M,K)
    y = y - y.min()
    yf2 = LA.norm(y,'fro')
    g = np.diag((np.transpose(D)@y)@Ek)
    
    
    lambdaval = np.max(g)/lambdadiv
    g = g - lambdaval
    
    #display data    
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    ax.plot(us)  # Plot some data on the axes.
    
    
    intSolver   = int(1)
    maxIter     = int(10)
    maxIterInt  = int(1000)
    addSize     = int(5)
        
    f1 = DSv7(yf2,D,Ek,lambdaval,H,g,intSolver,maxIter,maxIterInt,addSize)

    print("f1 = %f" %f1)
    return 

if __name__ == '__main__':
    print("Starting TB ....")
    DSv7TB()
    
