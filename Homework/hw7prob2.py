import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
import math
import matplotlib.pyplot as plt



def driver():

    N = 15
    h = 2 / (N-1)

    f = lambda x: 1 / (1 + (10*x)**2)
    xvals = [-1 + (j-1)* h for j in range(1,N+1)]

    yvals = [f(x) for x in xvals]


    
    ''' create points for evaluating the Lagrange interpolating polynomial'''
    Neval = 1000
    xeval = np.linspace(-1,1,Neval+1)
    yeval_l= np.zeros(Neval+1)


    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
        yeval_l[kk] = eval_lagrange(xeval[kk],xvals,yvals,N)



    ''' create vector with exact values'''
    fex = f(xeval)
    plt.figure()
    plt.plot(xeval,fex, color = 'steelblue', label = "True Function")
    plt.plot(xeval,yeval_l, color = 'red', label = "Barycentric Approximation")
    plt.plot(xvals, f(np.array(xvals)), 'o', label = "Interpolation Nodes")
    plt.title("Barycentric approximation using uniform nodes")

    plt.legend()
    plt.show()


def eval_lagrange(xeval,xint,yint,N):
    if xeval in xint:
        return yint[xint.index(xeval)]
    
    I = 1
    for n in range(N):
        I = I*(xeval - xint[n])
    mid_sum = 0
    for j in range(N):
        w = 1
        for i in range(N):
            if i != j:
                w = w/(xint[j]-xint[i])
        mid_sum += w * yint[j] / (xeval - xint[j])

    
    return I * mid_sum
    
driver()