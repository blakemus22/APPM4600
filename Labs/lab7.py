import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
import math
import matplotlib.pyplot as plt



def driver():

    N = 19
    h = 2 / (N-1)

    f = lambda x: 1 / (1 + 10*x**2)
    xvals = [-1 + (j-1)* h for j in range(N+1)]

    yvals = [f(x) for x in xvals]


    
    ''' create points for evaluating the Lagrange interpolating polynomial'''
    Neval = 1000
    xeval = np.linspace(-1,1,Neval+1)
    yeval_l= np.zeros(Neval+1)
    yeval_dd = np.zeros(Neval+1)
    '''Initialize and populate the first columns of the
    divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )
    for j in range(N+1):
        y[j][0] = yvals[j]

    y = dividedDiffTable(xvals, y, N+1)

    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
        yeval_l[kk] = eval_lagrange(xeval[kk],xvals,yvals,N)
        yeval_dd[kk] = evalDDpoly(xeval[kk],xvals,y,N)

    '''evaluate with monomial'''
    V = Vandermonde(xvals,N)
    Vinv = inv(V)
    coef = Vinv @ yvals
    yeval_mn = eval_monomial(xeval,coef,N,Neval)

    ''' create vector with exact values'''
    fex = f(xeval)
    plt.figure()
    plt.plot(xeval,fex, label = "True Function")
    plt.plot(xeval,yeval_l, label = "Lagrange")
    plt.plot(xeval,yeval_dd, label = "Newton DD")
    plt.plot(xeval, yeval_mn, label = "Monomial")
    plt.legend()


    plt.figure()
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    err_mn = abs(yeval_mn - fex)
    plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    plt.semilogy(xeval, err_dd,'bs--',label='Newton DD')
    plt.semilogy(xeval, err_mn, 'gv--', label = 'Monomial')
    plt.legend()
    plt.show()


def eval_lagrange(xeval,xint,yint,N):
    lj = np.ones(N+1)
    for count in range(N+1):
        for jj in range(N+1):
            if (jj != count):
                lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])
    yeval = 0.
    for jj in range(N+1):
        yeval = yeval + yint[jj]*lj[jj]
    return(yeval)


''' create divided difference matrix'''
def dividedDiffTable(x, y, n):
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /(x[j] - x[i + j]))
    return y


def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)
    ptmp[0] = 1.
    for j in range(N):
        ptmp[j+1] = ptmp[j]*(xval-xint[j])
    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N+1):
        yeval = yeval + y[0][j]*ptmp[j]
    return yeval

def eval_monomial(xeval,coef,N,Neval):
    yeval = coef[0]*np.ones(Neval+1)
    for j in range(1,N+1):
        for i in range(Neval+1):
            yeval[i] = yeval[i] + coef[j]*xeval[i]**j
    return yeval


def Vandermonde(xint,N):
    V = np.zeros((N+1,N+1))
    ''' fill the first column'''
    for j in range(N+1):
        V[j][0] = 1.0
    for i in range(1,N+1):
        for j in range(N+1):
            V[j][i] = xint[j]**i
    return V



driver()


# summary of the rest of the lab:
# we want to compare the three methods that we have looked at so far to create polynomial approximations 
# of functions. We can can plot the function in the driver on the interval [-1,1] with different degrees
# of polynomials. 
# what happens if we use different spacing of the intervals? how does this change the stability?

