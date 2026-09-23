# Working through hw2 as practice problems (was a DNF)
# Benjamin Ryan
import numpy as np
import matplotlib.pyplot as plt
A = np.array([[0.5,0.5],[(1+10**-10)/2, (1-10**-10)/2]])

A_inv = np.array([[1-10**10, 10**10],[1+10**10,-10**10]])


b = np.array([1,1])
x = np.array([1,1])

kA = np.linalg.norm(A)*np.linalg.norm(A_inv)
print('condition #:',kA)

x = 9.999999995e-10
res = []
nfac = 1
for n in range(1,10):
    nfac *= n
    res.append((np.e**x * x**(n+1))/(10**-9 * nfac))

print(res)

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

f = lambda x: x - 4*np.sin(2*x)-3

f_fixpt = lambda x: -1*np.sin(2*x) + 5*x/4 - 0.75


x = np.linspace(-2,8,100)

y = f(x)

[pts,ier] = fixedpt_iteration(f_fixpt,0,1e-10,50)
print(pts,ier)


plt.plot(x,y)
plt.hlines(0,-2,8,colors='black')
plt.show()

