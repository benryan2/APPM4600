import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math
from scipy.integrate import quad

def driver():

#  function you want to approximate
    f = lambda x: math.exp(x)

# Interval of interest    
    a = -1
    b = 1
# weight function    
    w = lambda x: 1.

# order of approximation
    n = 2

#  Number of points you want to sample in [a,b]
    N = 100
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])
      
    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])
        
    plt.figure()    
    plt.plot(xeval,fex, label= 'f(x)')
    plt.plot(xeval,pval,label= 'Expansion') 
    plt.legend()
    plt.show()    
    
    err = abs(pval-fex)
    # plt.semilogy(xeval,err_l,'ro--',label='error')
    plt.legend()
    plt.show()
    
def eval_legandre(n,x):
    p = [1,x]
    while len(p) != n+1:
        p.append( (1/(2*n + 1))*((2*n+1)*x*p[-1]-n*p[-2]))
    return p
    

def eval_legendre_expansion(f,a,b,w,n,x): 

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab 
  p = eval_legandre(n,x)
  # initialize the sum to 0 
  pval = 0.0    
  for j in range(0,n+1):
    # make a function handle for evaluating phi_j(x)
    phi_j = lambda x: eval_legandre(j,x)[j]
    # make a function handle for evaluating phi_j^2(x)*w(x)
    phi_j_sq = lambda x: (phi_j(x)**2) * w(x)
    # use the quad function from scipy to evaluate normalizations
    norm_fac,err = quad(phi_j_sq,a,b)
    # make a function handle for phi_j(x)*f(x)*w(x)/norm_fac
    func_j = lambda x: phi_j(x)*f(x)*w(x)
    # use the quad function from scipy to evaluate coeffs
    num,errn = quad(func_j,a,b)
    aj = norm_fac/num
    # accumulate into pval
    pval = pval+aj*phi_j(x)
       
  return pval
    
if __name__ == '__main__':
  # run the drivers only if this is called from the command line
  driver()               
