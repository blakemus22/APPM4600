import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
import math
import matplotlib.pyplot as plt



def driver():

    N = 17
    h = 2 / (N-1)

    f = lambda x: 1 / (1 + (10*x)**2)
    xint = [-1 + (j-1)* h for j in range(1,N+1)]

    yint = [f(x) for x in xint]

    Neval = 1000
    xeval = np.linspace(-1,1,Neval+1)

    '''evaluate with monomial'''
    V = Vandermonde(xint,N)
    Vinv = inv(V)
    coef = Vinv @ yint
    yeval_mn = eval_monomial(xeval,coef,N,Neval)

    ''' create vector with exact values'''
    fex = f(np.array(xint))
    plt.figure()
    plt.plot(xeval, f(xeval), color = 'steelblue', label = "True Function")
    plt.plot(xint,fex,'o', label = "Interpolating Nodes")
    plt.plot(xeval, yeval_mn, color = 'green', label = "Monomial Approximation")
    plt.title("Monomial Approximation using evenly spaced nodes")
    plt.legend()

    plt.show()








def eval_monomial(xeval,coef,N,Neval):
    yeval = coef[0]*np.ones(Neval+1)
    for j in range(1,N):
        for i in range(Neval+1):
            yeval[i] = yeval[i] + coef[j]*xeval[i]**j
    return yeval


def Vandermonde(xint,N):
    V = np.zeros((N,N))
    ''' fill the first column'''
    for j in range(N):
        V[j][0] = 1.0
    for i in range(1,N):
        for j in range(N):
            V[j][i] = xint[j]**i
    return V



driver()


