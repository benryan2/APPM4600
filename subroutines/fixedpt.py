# import libraries
import numpy as np
    
def driver():

# test functions 
     f1 = lambda x: 1+0.5*np.sin(x)
# fixed point is alpha1 = 1.4987....

     f2 = lambda x: 3+2*np.sin(x)
#fixed point is alpha2 = 3.09... 

     Nmax = 100
     tol = 1e-6

# test f1 '''
     x0 = 0.0
     [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
    
#test f2 '''
     x0 = 0.0
     [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f2(xstar):',f2(xstar))
     print('Error message reads:',ier)

#test f1 for iterations
     x0 = 0.0
     [pts,ier] = fixedpt_iteration(f1,x0,tol,Nmax)
     print('Vals at iterations',pts)
     print('f1(xstar):',f1(pts[-1]))
     print('Error message reads:',ier)

#test 2 for iterations
     x0 = 0.0
     [pts,ier] = fixedpt_iteration(f2,x0,tol,Nmax)
     print('Vals at iterations',pts)
     print('f1(xstar):',f1(pts[-1]))
     print('Error message reads:',ier)

# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    
def fixedpt_iteration(f,x0,tol,Nmax):
     '''
     x0 = initial guess
     tol = stopping tolerance
     nMax = max iterations
     '''
     points = np.zeros((Nmax,1))
     points[0] = x0
     for i in range(1,Nmax):
         points[i] = f(points[i-1])
         if(abs(points[i]-points[i-1]) < tol):
             points = points[:i+1]
             ier = 0
             return [points,ier]
         
     ier = 1
     return[points,ier]
    
# driver()
