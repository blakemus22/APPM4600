import numpy as np
import math
from numpy.linalg import norm


def driver():
    x0 = np.array([1,1,1])
    Nmax = 100
    tol = 1e-10
    [xstar,ier,its] = Surface(x0,tol,Nmax)
    print(xstar)
    print('Surface iteration: the error message reads:',ier)
    print('Surface iteration: number of iterations is:',its)


def evalF(x):
    F = np.zeros(3)
    denom = (2*x[0])**2 + (8*x[1])**2 + (8*x[2])**2
    f = x[0]**2 + 4*x[1]**2 + 4*x[2]**2 - 16
    F[0] = x[0] - (f * 2*x[0]) / denom
    F[1] = x[1] - (f * 8*x[1]) / denom
    F[2] = x[2] - (f * 8*x[2]) / denom
    return F



def Surface(x0,tol,Nmax):

    for its in range(Nmax):
        x1 = evalF(x0)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

driver()