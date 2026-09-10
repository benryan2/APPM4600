prelab 3

2.2 1)

Given a fixed point $p$ and the vector $\hat{p}$ of approximations by iteration, find a way to determine the order of convergence.

Considering the limit formula for convergence, a series converges with order $\alpha$ when 

$|p_{n+1}-p|\approx \lambda|p_n-p|^\alpha$

The terms $|p_{n+1}-p|$ and $|p_n-p|$ are errors between current iterations and the next iteration.

From a vector of approximations, we can subtract the actual value $p$ from each element to get an error vector $\hat{e}$.

Using the error vector, we then iterate over the vector to create an additional vector $\hat{\lambda}$, if this converges to a value between 0 and 1, the order of convergence is linear. You can check for convergence by establishing a tolerence.

If it does not converge, repeat, dividing the $nth$ term of $\hat{\lambda}$ by the $n+1th$ term of $\hat{e}$. check for convergence, repeat if necessary


the following code block is ugly.
```python
def convOrder(p_vector,tol,max_order=5):
    p = p_vector[-1]
    e = p_vector-p
    print(e)
    order = 0
    lam = e.copy()
    while(order < max_order):
        order += 1
        print(order)
        lam[0] /= e[1]
        for i in range(1, len(lam)-1):
            lam[i] = abs(lam[i] / e[i+1])
            if(abs(lam[i]-lam[i-1]) < tol):
                if(order == 1 and (lam[i] < 0 or lam[i] > 10)): break
                return [order,lam[i]]

    print('error')
    return [order,lam[i]]
```

2.2 2a)

It took 12 iterations to reach a fixed point with a tolerance of $10^{-10}$

2b)

Using the above method, I found that the algorithm converges at a 2nd order rate, and the associated value of $\lambda$ is $\lambda\approx7.85952371$