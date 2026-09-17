**Benjamin Ryan Lab 4 APPM4600**

Exercises: 3.1

Derive a condition which guaruntees that Newtons method will converge to a unique root for all initial guess in a neighborhood of the root.

This occurs when $|g'(x_0)| < 1$

where $g(x) = x - \frac{f(x)}{f'(x)}$

so $g'(x) = -\frac{f(x)f''(x)}{(f'(x))^2}$

This is implemented as a check before performing newton's method
```python
    gpx = abs(f(p0)*fpp(p0)/(fp(p0)**2))
    if gpx < 1:
        newton(f, fp, p0,tol,Nmax)
```

2) Bisection will cut between a and b, approaching the root until the midpoint is in the basin of convergence

3) the inputs to the bisection needed to include f, fp, and fpp. that is, the function and two of it's derivatives. It also needed the tolerance and Nmax values to be sent into newton's method

4) This method was combined as the function 'bisect_newton'

```python
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
```

5) The advantage of the new method is that it guaruntees that newton's method will converge, it's less likely to fail, like the bisection method, but faster, since it switches to newton's method once possible, accelerating to 2nd order convergence. It's limitation is that now you must have an interval in which there is a root and a sign change. The input function must also have two continuous derivatives, which must be provided.

