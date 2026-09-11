import numpy as np
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from subroutines.fixedpt import fixedpt_iteration


def driver():
    g = lambda x: (10/(x+4))**0.5
    p0 = 1.5
    tol = 10**-10
    nmax = 1000
    [pts, ier] = fixedpt_iteration(g,p0,tol,nmax)
    print(pts,ier)
    if ier == 1: return

    [order,lam] = convOrder(pts,10**-3)
    print(order,lam)

    alpha = convOrder_analytical(pts)
    print('analytically found alpha:',alpha)

def convOrder(p_vector,tol,max_order=5):
    p = p_vector[-1]
    e = abs(p_vector-p)
    print(e)
    order = 0
    lam = e.copy()
    while(order < max_order):
        order += 1
        print(order)
        lam[0] /= e[1]
        for i in range(1, len(lam)-1):
            lam[i] = abs(lam[i] / e[i+1])
            print(lam[i])
            if(abs(lam[i]-lam[i-1]) < tol):
                if(order == 1 and (lam[i] < 0 or lam[i] > 1)): break
                return [order,lam[i]]

    print('error')
    return [order,lam[i]]

def convOrder_analytical(p_vector):
    pstar = p_vector[-1]
    e_k = abs(p_vector[-2] - pstar)
    e_k1 = abs(p_vector[-3] - pstar)

    alpha = abs(np.log(e_k1/e_k)/np.log(e_k/e_k1))
    return alpha

# driver()