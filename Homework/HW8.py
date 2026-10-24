# Benjamin Ryan
# APPM 4600 HW8

import numpy as np
import matplotlib.pyplot as plt

def legrange_poly(pts,j,x):
    # returns legrange polynomial j evaluated at x (vector)
    result = np.ones(len(x))
    for i in range(len(pts[0])):
        if i == j: continue
        result *= ((x - pts[0][i])/(pts[0][j] - pts[0][i]))
    return result

def legrange_derivative(pts,j,x):
    # returns the first derivative of the legrange polynomial j eval at x
    sigma = 0
    for i in range(len(pts[0])):
        if i == j: continue
        sigma += (x-pts[0][i])**-1
    return legrange_poly(pts,j,x)*sigma


def legrange_interp(pts,n_x):
    # returns x,y
    x = np.linspace(pts[0][0],pts[0][-1],n_x)
    y = np.zeros(len(x))
    for j in range(len(pts[0])):
        y += (pts[1][j]*legrange_poly(pts,j,x))

    return x,y

def hermite(pts,n_x):
    x = np.linspace(pts[0][0],pts[0][-1],n_x)  
    y = np.zeros(len(x))
    for j in range(len(pts[0])):
        l = legrange_poly(pts,j,x)
        lp = legrange_derivative(pts,j,x)
        Q = (1 - 2*(x-pts[0][j])*lp)*(l**2)
        R = (x-pts[0][j])*(l**2)
        y += (pts[1][j]*Q + pts[2][j]*R)
    return x,y

def natural_cubic(pts,n_x):
    n = len(pts[0])
    order = np.argsort(pts[0])
    x_p = pts[0][order]
    y_p = pts[1][order]

    # Calculate Mk values
    Y = np.zeros(n)
    A = np.zeros((n,n))
    A[0][0]=1
    A[-1][-1]=1
    for k in range(1,n-1):
        Y[k] = ((y_p[k+1]-y_p[k])/(x_p[k+1]-x_p[k])
                - (y_p[k]-y_p[k-1])/(x_p[k]-x_p[k-1]))
        A[k][k-1] = (x_p[k]-x_p[k-1])/6
        A[k][k+1] = (x_p[k+1]-x_p[k])/6
        A[k][k] = ((x_p[k]-x_p[k-1])+(x_p[k+1]-x_p[k]))/3

    M = np.linalg.inv(A) @ Y
    print(Y)
    print(A)
    print(M)

    # Calculate Ck, Dk, spline
    x = np.linspace(x_p[0],x_p[-1],n_x)
    y = np.zeros(len(x))
    
    for k in range(n-1):
        hk = x_p[k+1] - x_p[k]
        Ck = y_p[k]/hk - (M[k]*hk)/6
        Dk = y_p[k+1]/hk - (M[k+1]*hk)/6
        section = (x_p[k] <= x) & (x <= x_p[k+1])
        y[section] = ( (M[k]/(6*hk))*(x_p[k+1]-x[section])**3
                      + (M[k+1]/(6*hk))*(x[section]-x_p[k])**3
                      + Ck*(x_p[k+1]-x[section]) + Dk*(x[section]-x_p[k]) )

    return x,y

def clamped_cubic(pts,n_x):
    n = len(pts[0])
    order = np.argsort(pts[0])
    x_p = pts[0][order]
    y_p = pts[1][order]

    # Calculate Mk values
    Y = np.zeros(n)
    A = np.zeros((n,n))
    A[0][0]= (x_p[1] - x_p[0])/3
    A[0][1] = A[0][0]/2
    A[-1][-1]= (x_p[-1]-x_p[-2])/3
    A[-1][-2] = A[-1][-1]/2
    for k in range(1,n-1):
        Y[k] = ((y_p[k+1]-y_p[k])/(x_p[k+1]-x_p[k])
                - (y_p[k]-y_p[k-1])/(x_p[k]-x_p[k-1]))
        A[k][k-1] = (x_p[k]-x_p[k-1])/6
        A[k][k+1] = (x_p[k+1]-x_p[k])/6
        A[k][k] = ((x_p[k]-x_p[k-1])+(x_p[k+1]-x_p[k]))/3

    M = np.linalg.inv(A) @ Y
    print(Y)
    print(A)
    print(M)

    # Calculate Ck, Dk, spline
    x = np.linspace(x_p[0],x_p[-1],n_x)
    y = np.zeros(len(x))
    
    for k in range(n-1):
        hk = x_p[k+1] - x_p[k]
        Ck = y_p[k]/hk - (M[k]*hk)/6
        Dk = y_p[k+1]/hk - (M[k+1]*hk)/6
        section = (x_p[k] <= x) & (x <= x_p[k+1])
        y[section] = ( (M[k]/(6*hk))*(x_p[k+1]-x[section])**3
                        + (M[k+1]/(6*hk))*(x[section]-x_p[k])**3
                        + Ck*(x_p[k+1]-x[section]) + Dk*(x[section]-x_p[k]) )

    return x,y

def periodic_cubic(pts,n_x):
    n = len(pts[0])
    order = np.argsort(pts[0])
    x_p = pts[0][order]
    y_p = pts[1][order]

    # Calculate Mk values
    Y = np.zeros(n)
    A = np.zeros((n,n))
    h0 = x_p[1]-x_p[0]
    hn_1 = x_p[-1] - x_p[-2]
    A[0][0] = (h0+hn_1)/3
    A[0][1] = h0/6
    A[0][-2] = hn_1/6

    A[-1][-1] = (h0+hn_1)/3
    A[-1][1] = h0/6
    A[-1][-2] = hn_1/6

    Y[0] = (y_p[1]-y_p[0])/h0 - (y_p[-1]-y_p[-2])/hn_1
    Y[-1] = (y_p[1]-y_p[0])/h0 - (y_p[-1]-y_p[-2])/hn_1

    for k in range(1,n-1):
        Y[k] = ((y_p[k+1]-y_p[k])/(x_p[k+1]-x_p[k])
                - (y_p[k]-y_p[k-1])/(x_p[k]-x_p[k-1]))
        A[k][k-1] = (x_p[k]-x_p[k-1])/6
        A[k][k+1] = (x_p[k+1]-x_p[k])/6
        A[k][k] = ((x_p[k]-x_p[k-1])+(x_p[k+1]-x_p[k]))/3

    M = np.linalg.inv(A) @ Y
    print(Y)
    print(A)
    print(M)

    # Calculate Ck, Dk, spline
    x = np.linspace(x_p[0],x_p[-1],n_x)
    y = np.zeros(len(x))
    
    for k in range(n-1):
        hk = x_p[k+1] - x_p[k]
        Ck = y_p[k]/hk - (M[k]*hk)/6
        Dk = y_p[k+1]/hk - (M[k+1]*hk)/6
        section = (x_p[k] <= x) & (x <= x_p[k+1])
        y[section] = ( (M[k]/(6*hk))*(x_p[k+1]-x[section])**3
                      + (M[k+1]/(6*hk))*(x[section]-x_p[k])**3
                      + Ck*(x_p[k+1]-x[section]) + Dk*(x[section]-x_p[k]) )

    return x,y


def get_nodes(f,fp,a,b,n):
    x = np.linspace(a,b,n)
    y = f(x)
    yp = fp(x)
    return x,y,yp

def get_chebyshev(f,fp,a,b,n):
    x_i = []
    for i in range(n):
        i = i+1
        x_i.append(np.cos((np.pi*(2*i - 1))/(2*n)))
    x = np.array(x_i)
    x[x<0]*=abs(a)
    x[x>0]*=abs(b)
    y = f(x)
    yp = fp(x)
    return x,y,yp


def driver():
    f = lambda x: 1/(1+x**2)
    fp = lambda x: (-2*x)/((1+x**2)**2)
    n_vals = [5,10,15,20]

    '''
    # Legrange
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,-5,5,n)
        x,y = legrange_interp(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Legrange Interpolation From Evenly Spaced Points")
    plt.show() 

    # Hermite
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,-5,5,n)
        x,y = hermite(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Hermite Interpolation From Evenly Spaced Points")
    plt.show()

    # cubic spline
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,-5,5,n)
        x,y = natural_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Natural Cubic Spline From Evenly Spaced Points")
    plt.show()
    
    # clamped cubic spline
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,-5,5,n)
        x,y = clamped_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Clamped Cubic Spline From Evenly Spaced Points")
    plt.show()

    # Legrange Cheby
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_chebyshev(f,fp,-5,5,n)
        x,y = legrange_interp(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Legrange Interpolation From Chebyshev Points")
    plt.show() 

    # Hermite cheb
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_chebyshev(f,fp,-5,5,n)
        x,y = hermite(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Hermite Interpolation From Chebyshev Points")
    plt.show()

    # cubic spline cheb
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_chebyshev(f,fp,-5,5,n)
        x,y = natural_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Natural Cubic Spline From Chebyshev Points")
    plt.show()
    
    # clamped cubic spline chev
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_chebyshev(f,fp,-5,5,n)
        x,y = clamped_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Clamped Cubic Spline From Chebyshev Points")
    plt.show()
    '''
    # Periodic Spline
    f = lambda x: np.sin(10*x)
    fp = lambda x: 10*np.cos(10*x)

    ig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,0,2*np.pi,n)
        x,y = natural_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Natural Cubic Spline From Evenly Spaced Points")
    plt.show()

    
    ig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,8))
    for i,n in enumerate(n_vals):
        pts = get_nodes(f,fp,0,2*np.pi,n)
        x,y = periodic_cubic(pts,1001)
        ax[i//2,i%2].scatter(pts[0],pts[1],color='orange')
        ax[i//2,i%2].plot(x,y)
        ax[i//2,i%2].set_title(f'{n} Points')

    plt.suptitle("Periodic Spline From Evenly Spaced Points")
    plt.show()


driver()