# Benjamin Ryan
'''
Lab 4: Creating more robust root finding methods


'''
import numpy as np

def driver():
    f = lambda x: np.e**(x**2 + 7*x - 30) - 1
    fp = lambda x: (2*x + 7)*np.e**(x**2 + 7*x - 30)
    fpp = lambda x: (4*x**2 + 14*x + 51)*np.e**(x**2 + 7*x - 30)
    # expected root = 3
    tol = 10**-5
    Nmax = 100

    [newtonP,newtonPstar,ier,m] = newton(f, fp, 4.5, tol, Nmax)
    


def bisect_newton(f,fp,fpp,a, b,tol,Nmax):
    #confirm f(a) and f(b) have opposite signs
    if f(a)*f(b) >= 0:
        ier = 1
        return [[], a, ier]

    n = 0
    while n < Nmax:
        mid = 0.5*(a+b)
        n += 1
        #check if mid is root
        if f(mid) == 0:
            ier = 0
            return [[],mid,ier]
        # check if mid in basin
        gpx = abs(f(mid)*fpp(mid)/(fp(mid)**2))
        if gpx < 1:
            return newton(f, fp, mid, tol, Nmax);
        #if not in mid, iterate
        if f(a)*f(mid) < 0:
            b = mid
        else:
            a = mid


def newton(f, fp, p0, tol, Nmax):
    p = np.zeros(Nmax+1);
    p[0] = p0
    for it in range(Nmax):
        p1 = p0-f(p0)/fp(p0)
        p[it+1] = p1
        if (abs(p1-p0) < tol):
            pstar = p1
            info = 0
            return [p[:it+1],pstar,info,it]
        p0 = p1
    pstar = p1
    info = 1
    return [p,pstar,info,it]
        
driver()
