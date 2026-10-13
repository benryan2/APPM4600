# Benjamin Ryan 
# APPM 4600 Homework 7
import numpy as np
import matplotlib.pyplot as plt

def coeff_interp(f,x_pts,N,x):
    p = np.arange(N)
    y_pts = f(x_pts)
    
    X = np.repeat(np.array([x_pts]),N,axis=0).T
    P = np.repeat(np.array([p]),N,axis=0)
    V = X ** P

    print('X\n',X)
    print('P\n',P)
    print('V\n',V)

    c = np.linalg.inv(V) @ y_pts
    return (np.repeat([x],N,axis=0).T**np.repeat([np.arange(N)],len(x),axis=0))@c
    # note to self: wtf is this ^ just use a loop, asshole.
        # thank god that works, i wrote that high as hell

def phi(x,x_vec,N):
    val = 1
    for i in range(N):
        val *= (x - x_vec[i])
    return val

def wj(x_vec,j,N):
    val = 1
    for i in range(N):
        if i == j: continue
        val *= (x_vec[j]-x_vec[i])
    return 1/val
    
def bary_leg(f,x_vec,xrange):
    # returns a vector y = p(xrange)
    n = len(x_vec)
    y = []
    for index,x in enumerate(xrange):
        sum_val = 0
        for j in range(n):
            sum_val += wj(x_vec,j,n)/(x - x_vec[j])
        y.append(phi(x,x_vec,n)*sum_val*f(x))

    return np.array(y)
    
def problem1(N):
    f = lambda x: 1/(1+(10*x)**2)
    x_pts = np.linspace(-1,1,N)
    print(x_pts)
    print(f(x_pts))

    
    x = np.linspace(-1,1,1001)
    y = coeff_interp(f,x_pts,N,x)
    
    print(max(y))

    # Plot
    plt.plot(x,y)
    plt.scatter(x_pts,f(x_pts),marker='o',color='orange')
    # plt.plot(x,f(x))
    plt.title(f"Interpolated Polynomial From {N} Points")
    plt.show()

def problem2(N):
    # Repeat interpolation of f(x) using barycentric Lagrange interpolation
    f = lambda x: 1/(1 + (10*x)**2)
    x_vec = np.linspace(-1,1,N)
    print(x_vec)

    x = np.linspace(-1,1,1001)
    actual_y = f(x)
    y = bary_leg(f,x_vec,x)
    err = abs(y - actual_y)

    plt.scatter(x_vec,f(x_vec),color='orange')
    plt.plot(x,y)
    plt.title(f"Barycentric Legrange Interpolation from {N} Points")
    plt.show()

    plt.semilogy(x,err)
    plt.title(f"Absolute Error in Legrange Interpolation from {N} Points")
    plt.show()

def problem3(N):
    f = lambda x: 1/(1+(10*x)**2)
    x_i = []
    for i in range(N):
        i = i+1
        x_i.append(np.cos((np.pi*(2*i - 1))/(2*N)))
    x_vec = np.array(x_i)

    print(x_vec)

    x = np.linspace(-1,1,1001)
    y = f(x)
    y_coef = coeff_interp(f,x_vec,N,x)
    y_bary = bary_leg(f,x_vec,x)
    err_bary = abs(y - y_bary)

    plt.plot(x_vec, f(x_vec),'o')
    # plt.plot(x,y_coef)
    plt.plot(x,y_bary)
    plt.title(f"Legrange Interpolation From {N} Chebyshev Points")
    plt.show()

    plt.semilogy(x,err_bary)
    plt.title(f"Absolute Error in Legrange Interpolation From {N} Chebyshev Points")
    plt.show()

    plt.plot(x_vec,f(x_vec),'o')
    plt.plot(x,y_coef)
    plt.title(f"Coefficient Interpolation From {N} Chebyshev Points")
    plt.show()


# problem1(19)
# problem2(150)
problem3(101)