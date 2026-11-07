# Benjamin Ryan
# APPM 4600
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def problem1():
    I = -12.070346316389632
    f = lambda x: (np.e**x)*np.cos(x)
    a = 0; b = np.pi
    tol = 1e-10
    trap_val, h_trap, err_trap = trap_to_tol(f,a,b,I,tol,100)
    simp_val, h_simp, err_simp = simpson_to_tol(f,a,b,I,tol,100)

    plt.loglog(h_trap,err_trap)
    plt.loglog(h_simp,err_simp)
    plt.xlabel('Subinterval Length h')
    plt.ylabel('Absolute Error')
    plt.title('Convergence of Trapezoid and Simpsons Quadrature Methods')
    plt.show()




def trap_to_tol(f,a,b,I,tol,max_it):
    n = 2
    it=0
    A = np.inf
    err = []
    h_list = []
    while it < max_it:
        it += 1
        x = np.linspace(a,b,n)
        area = 0
        h = (b-a)/(n-1)
        for i in range(len(x)-1):
            area += 0.5*h*(f(x[i])+f(x[i+1]))
        print(it,n,h,area)
        err.append(abs(area-I))
        h_list.append(h)
        if abs(area-I) < tol:
            return area, h_list, err
        else:
            n = (n*2)-1
    print('err max it')

def simpson_to_tol(f,a,b,I,tol,max_it):
    n = 3
    it = 0 
    A = np.inf
    err = []
    h_list = []
    while it < max_it:
        it += 1
        x = np.linspace(a,b,n)
        area = 0
        h = (b-a)/(n-1)
        i = 0
        while i < len(x)-2:
            area += (h/3)*(f(x[i])+4*f(x[i+1])+f(x[i+2]))
            i+=2
        print(it,n,h,area)
        err.append(abs(area-I))
        h_list.append(h)
        if abs(area-I) < tol:
            return area, h_list, err
        else:
            n = (n*2)-1
        
def trap_set_n(f,a,b,n):
    x = np.linspace(a,b,n)
    I = 0
    h = (b-a)/(n-1)
    for i in range(len(x)-1):
        I += (h/2)*(f(x[i])+f(x[i+1]))

    return I


def problem2():
    f = lambda s: 1/(1+s**2)
    a = -5; b = 5

    quadres1, err1 = quad(f,a,b)
    quadres2, err2, infodict2 = quad(f,a,b,epsabs=1e-4,full_output=True)

    T_n = trap_set_n(f,a,b,190)
    quad
    ttol1 = abs(T_n-quad(f,a,b)[0])
    ttol2 = abs(T_n-quad(f,a,b,epsabs=1e-4)[0])
    print("\nTRAPEZOID:")
    print(f"T_n for n=190: {T_n}")
    print(f"Absolute error to Scipy quad, default tolerance: {ttol1}")
    print(f"Absolute error to Scipy quad, tolerance=10^-4: {ttol2}")
    print(f"Number of function evals by Scipy Quad for 10^-4 tolerance: {infodict2['neval']}")











# problem1()

problem2()