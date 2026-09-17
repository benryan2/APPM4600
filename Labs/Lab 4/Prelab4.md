Benjamin Ryan

Lab 4 prelab

Comparing methods for root finding

- Bisection
    - 
    - Input
        - f, a, b, tolerance
    - Iteration
        - checks midpoint, restarts on new interval based on sign of midpoint
    - Idea behind method
        - binary search
    - Required for convergence
        - there must be a sign change between a and b,
        - funciton must be continuous on (a,b)
    - Pros
        - fast
    - Cons 
        - Requires knowing more or less where root is, needs a sign change

- Fixed point
    -
    - Input
        - f, a, b, tolerance, max_i
    - Iteration
        - $x_{n+1} = f(x_n)$
    - Idea behind method
        - rapid, simple iterations to finding where $f(x) = x$
        - find root by making a new funciton $g(x) = f(x) + x$
    
    - Required for Convergence
        - $f$ must be continuous
        - $f(x)$ must be between a,b
        - $f'(x)$ must exist
    - pros
        - Simple, reliable
    - cons
        - if $f(x) < b$ at any point, fails to converge

- Newton
    -
    - Input
        - $a, b, x_0,$ max_iter, tol
    - Iteration
        - $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
    - Idea behind method
        - Follows tangent lines to x-intercept, eventually will be a root
    - Required for Convergence
        - $f'(x)$ must exist
        - $f''(x)$ must exist
        - $f(x_n) \neq 0$ for any $x_n$ 
    - pros
        - fast convergence
    - cons
        - can fail, needs 2 derivatives

- Secant
    -
    - Input
        - $x_0, x_1,$ tol, max_iter
    - Iteration
        - $x_{n+1} = x_n - f(x_n)\frac{(x_n - x_{n-1})}{f(x_n) - f(x_{n-1})}$
    - Idea behind method
        - uses the secant line to find x-intersection, similar to newton's method
    - Required for Convergence
        - secant line can't be horizontal
    - pros
        - does not require calculating derivative of $f(x)$
    - cons
        - Requires two points for input.
