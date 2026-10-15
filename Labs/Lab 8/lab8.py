# Benjamin Ryan Lab 8
# APPM 4600
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.linalg import inv 

def prelab(x0,y0,x1,y1,a):
    b = (y1-y0)/(x1-x0)
    return b*(a - x0) + y0

def driver():

    f = lambda x: 1/(1 + (10*x)**2)
    a = -1
    b = 1

    ''' create points you want to evaluate at '''
    Neval = 100
    xeval = np.linspace(a,b,Neval)

    ''' number of intervals '''
    Nint = 10

    '''evaluate the linear spline'''
    yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)

    ''' evaluate f at the evaluation points'''
    fex = f(xeval)

    plt.plot(xeval,fex,label='exact function')
    plt.plot(xeval,yeval,label='natural spline')
    plt.legend()
    plt.show()

    err = abs(yeval-fex)
    plt.figure()
    plt.semilogy(xeval,err,'ro--',label='absolute error')
    plt.legend()
    plt.show()

def  eval_lin_spline(xeval,Neval,a,b,f,Nint):

    '''create the intervals for piecewise approximations'''
    xint = np.linspace(a,b,Nint+1)
   
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval) 
    
    for jint in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        '''let n denote the length of ind'''
        
        '''temporarily store your info for creating a line in the interval of 
         interest'''
        a1= xint[jint]
        fa1 = f(a1)
        b1 = xint[jint+1]
        fb1 = f(b1)
        print(a1,b1)
        n = Neval//Nint

        for kk in range(n):
            '''use your line evaluator to evaluate the lines at each of the points 
            in the interval'''
            if a1 <= xeval[kk] <= b1:
                yeval[kk] = prelab(a1,fa1,b1,fb1,xeval[kk])
            '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
            the points (a1,fa1) and (b1,fb1)'''
            

    return yeval

driver()