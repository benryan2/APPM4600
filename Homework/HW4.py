# Benjamin Ryan
# APPM 4600 HW 4 work

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
          return [p[:it+2],pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def newton_m(f,fp,m,p0,tol,Nmax):
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
      p1 = p0- m* f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p[:it+2],pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def problem3():
    f = lambda x: (np.e**x - 3*x**2)**3
    fp = lambda x: 3*(np.e**x - 6*x)*(np.e**x - 3*x**2)**2
    g = lambda x: (np.e**x - 3*x**2)/(3*(np.e**x - 6*x))
    gp = lambda x: (1 - (np.e**x - 3*x**2)*(np.e**x - 6)*(np.e**x - 6*x)**-2) / 3

    # plot function
    x = np.linspace(3,5,100)
    y = f(x)
    plt.plot(x,y)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.show()

    [p_n,pstar_n,info_n,it_n] = newton(f,fp,4,10**-13,100)
    [p_2,pstar_2,info_2,it_2] = newton_m(f,fp,3,4,10**-13,100)
    [p_g,pstar_g,info_g,it_g] = newton(g,gp,4,10**-13,100)
    print(p_n,it_n,'\n',p_2,it_2,'\n',p_g,it_g)
    plt.semilogy(np.arange(len(p_n)-1),abs(p_n[:-1]-pstar_n),label='Newtons Method')
    plt.semilogy(np.arange(len(p_2)-1),abs(p_2[:-1]-pstar_2),label='Modified with m')
    plt.semilogy(np.arange(len(p_g)-1),abs(p_g[:-1]-pstar_g),label='Modified with g(x)',linestyle='--')

    plt.title('Error in Newtons Method Variations')
    plt.xlabel('Iteration n')
    plt.ylabel('$|p_n - p^*|$')
    plt.legend()
    plt.show()

# problem3()

def secant(f,p0,p1,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f - function
    p0,p1   - initial guesses
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+2);
  p[0] = p0
  p[1] = p1
  for it in range(Nmax):
      p2 = p1-f(p1)*(p1-p0)/(f(p1)-f(p0))
      p[it+2] = p2
      if (abs(p2-p1) < tol):
          pstar = p2
          info = 0
          return [p[:it+2],pstar,info,it]
      p0 = p1
      p1 = p2
  pstar = p2
  info = 1
  return [p,pstar,info,it]

def problem4():
    f = lambda x: x**6 - x - 1
    fp = lambda x: 6*x**5 - 1

    tol = 10**-16
    nmax = 100
    x0 = 2
    x1 = 1

    
    [p,pstar,ier,iter] = newton(f,fp,x0,tol,nmax)
    [s,sstar,sier,siter] = secant(f,x0,x1,tol,nmax)
    print(p,pstar,ier,iter)
    print(s,sstar,sier,siter)

    errp = p-pstar
    errs = s-sstar
    print(errp)
    print(errs)

    nrows = max(len(errp), len(errs))
    df = pd.DataFrame({
      'Newton Error': np.pad(errp, (0, nrows - len(errp)), constant_values=np.nan),
      'Secant Error': np.pad(errs, (0, nrows - len(errs)), constant_values=np.nan),
    })

    _, ax = plt.subplots(figsize=(8, 3))
    ax.axis('off')

    # Automatically converts the DataFrame into a Matplotlib table
    table = ax.table(
      cellText=df.values,
      colLabels=df.columns,
      rowLabels=df.index,
      loc='center',
      colWidths=[0.2, 0.2],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    table.scale(1, 1.5)

    plt.show()

    _,ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.plot(abs(errp[:-1]), abs(errp[1:]), label='Newton')
    plt.plot(abs(errs[:-1]), abs(errs[1:]), label='Secant')

    plt.title('Error Analysis for Newton and Secant Methods')
    plt.xlabel(r'$|x_k-\alpha|$')
    plt.ylabel(r"$|x_{k+1}-\alpha|$")
    plt.legend()
    plt.show()

    p_newton = (np.log(abs(errp[-3]/errp[-4])))/(np.log(abs(errp[-4]/errp[-5])))
    p_secant = (np.log(abs(errs[-2]/errs[-3])))/(np.log(abs(errs[-3]/errs[-4])))
    print("p Newton:",p_newton)
    print("p Secant:",p_secant)


problem4()