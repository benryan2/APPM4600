# Benjamin Ryan lab 3
# Numerical Methods APPM 4600

import numpy as np
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from subroutines.fixedpt import fixedpt_iteration
import find_order


def driver():
    g = lambda x: (10/(x+4))**0.5
    p0 = 1.5
    tol = 10**-10
    nmax = 1000
    [pts, ier] = fixedpt_iteration(g,p0,tol,nmax)
    if ier == 1: 
        print('fixedpoint failed')
        return

    print('points from fixedpt iterations:\n',pts)
    print('num iterations:',len(pts))

    p_hat = aitken(pts)
    print('Aitken output\n', p_hat)
    ord_ait = find_order.convOrder_analytical(p_hat)
    print('order = ',ord_ait)


def aitken(p):
    p_hat = []
    for i in range(len(p)-2):
        p_hat.append( p[i] - (p[i+1]-p[i])**2 / (p[i+2]-2*p[i+1]+p[i]))

    return np.array(p_hat)


driver()