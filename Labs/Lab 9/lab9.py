# Benjamin ryan APPM4600 lab 9
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def eval_legandre(n,x):
    p = [1,x]
    while len(p) != n:
        p.append( (1/(2*n + 1))*((2*n+1)*x*p[-1]-n*p[-2]))
    return p

def eval_fphiw(f,w,x,phi):
    return f(x)*phi[j]


def driver():
    quad()


driver()