Lab 3, Benjamin Ryan

3.1 Aitken's $\delta^2$ acceleration technique

Given a vector of values converging linearly to p, create a new sequence of approx, following

$\hat{p_n} = p_n - \frac{(p_{n+1}-p_n)^2}{p_{n+2}-2p_{n+1}+p_n}$

3.2 Exercises

```python
def aitken(p):
    p_hat = []
    for i in range(len(p)-2):
        p_hat.append( p[i] - (p[i+1]-p[i])**2 / (p[i+2]-2*p[i+1]+p[i]))

    return np.array(p_hat)
```
The convergence is in fact faster, reaching a steady state within 4 iterations.
