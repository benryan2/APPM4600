**Benjamin Ryan Lab 8 APPM4600**

***Prelab***

Write a subroutine that constructs a line that goes through $(x_0,f(x_0))$ and $(x_1,f(x_1))$, and evaluates at point $\alpha$

I completed this by calculating the slope, then solving the point slope formula.

```python
def prelab(x0,y0,x1,y1,a):
    b = (y1-y0)/(x1-x0)
    return b*(a - x0) + y0
```

***Lab Day***

Worked through 3, lots of code debuggin, code all on github!