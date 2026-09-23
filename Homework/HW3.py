#Benjamin Ryan HW3
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

def problem1():
    #problem 1
    fa = lambda x: -16 + 6*x + 12/x;
    fb = lambda x: (2/3)*x + x**-2
    fc = lambda x: 12/(1 + x)

    tol = 10**-8
    nMax = 100

    pa = evaluate_series(fa,1.8,tol,nMax)
    pb = evaluate_series(fb,1,tol,nMax)
    pc = evaluate_series(fc,1,tol,nMax)

    [alphaA,lamA] = convOrder_analytical(pa)
    [alphaB, lamB] = convOrder_analytical(pb)
    [alphaC, lamC] = convOrder_analytical(pc)

    print("series A:\n",pa,'convergence order:',alphaA,lamA,'\n\n')
    print("series B:\n",pb,'convergence order:',alphaB,lamB,'\n\n')
    print("series C:\n",pc,'convergence order:',alphaC,lamC,'\n\n')

def evaluate_series(f,x0,tol,nMax):
    p = [x0]
    p.append(f(x0))
    print(p)
    n = 1
    while n < nMax:
        p.append(f(p[-1]))
        n += 1
        if(abs(p[-1]-p[-2]) < tol): return p

    print('no convergence')
    return p

def convOrder_analytical(p_vector):
    for i in range(1,len(p_vector)-2):
        pstar = p_vector[-i]
        e_k = abs(p_vector[-i-1] - pstar)
        e_k1 = abs(p_vector[-i-2] - pstar)
        if e_k*e_k1 != 0: break

    alpha = abs(np.log(e_k1/e_k)/np.log(e_k/e_k1))
    lam = e_k1/e_k
    return [alpha,lam]

# problem1()

def problem2():
    Ti = 20
    Ts = -15
    alpha = 0.138e-6
    t = 40*24*60*60
    erf_integrand = lambda s: np.exp(-1*(s**2))
    
    def f(x):
        return Ts + (Ti - Ts) * erf(x / (2 * np.sqrt(alpha * t)))

    def df(x):
        Ti = 20
        Ts = -15
        alpha = 0.138e-6
        t = 60 * 24 * 60 * 60 
        
        return ((Ti - Ts) / np.sqrt(np.pi * alpha * t)) * np.exp(-(x**2) / (4 * alpha * t))

    x_vals = np.linspace(0, 3, 200) 
    y_vals = f(x_vals)
    
    plt.plot(x_vals, y_vals, label="f(x) = Temperature")
    plt.xlabel("Depth (meters)")
    plt.ylabel("Temperature (°C)")
    plt.title("Temperature at depth")
    plt.legend()
    plt.grid(True)
    plt.show()

    [root_bisection, err_b] = bisection(f,0,3,10**-13)
    print("root from bisection:",root_bisection)
    [p,pstar,info,it] = newton(f,df,0.01,10**-13,100)
    print("root from newton:",p)




def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    print('count = ', count)
    return [astar, ier]
                   
def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]


# problem2()

def problem3():
   fa = lambda x: x*(1+(7-x**5)/x**2)**3
   fb = lambda x: x - (x**5 - 7)/x**2
   fc = lambda x: x - (x**5 - 7)/(5*x**4)
   fd = lambda x: x - (x**5 - 7)/12

   [pa,era] = fixedpt_iteration(fa,1,10**-10,100)
   [pb,erb] = fixedpt_iteration(fb,1,10**-10,100)
   [pc,erc] = fixedpt_iteration(fc,1,10**-10,100)
   [pd,erd] = fixedpt_iteration(fd,1,10**-10,100)

   print("A:",pa,era,'\n\n')
   print("B:",pb,erb,'\n\n')
   print("C:",pc,erc,'\n\n')
   print("D:",pd,erd,'\n')

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

# problem3()
f = lambda x: x**3 + x - 4
[res,ier] = bisection(f,1,4,1e-3)
print(res,ier)