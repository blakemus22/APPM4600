import math
import numpy as np
import matplotlib.pyplot as plt


def driver():


    x0 = 4

    tol = 1e-13

    # problem 4
    f4 = lambda x: np.e**(3*x) -27*x**6 + 27*x**4*np.e**x - 9*x**2*np.e**(2*x)
    f4p = lambda x: 3*np.e**(3*x) -162*x**5 +27*x**4*np.e**x + 108*x**3*np.e**x-18*x**2*np.e**(2*x) -18*x*np.e**(2*x)
    f4pp = lambda x: 9*np.e**(3*x) -810*x**4 +27*x**4*np.e**x +216*x**3*np.e**x + 324*x**2*np.e**x -36*x**2*np.e**(2*x) -72*x*np.e**(2*x) -18*np.e**(2*x)


    [a,astar, ier, it] = newton(f4,f4p,f4pp, x0, tol, 300)
    print('the approximate root is', astar)
    print('the error message reads:', ier)
    print('f(astar)=', f4(astar))
    print(len(a), 'iterations to converge')


# define routines

def newton(f,fp, fpp, p0,tol,Nmax):

    p = np.zeros(Nmax+1)
    p[0] = p0
    for it in range(Nmax):
        # p1 = p0- f(p0)/fp(p0)
        # modified part ii
        # p1 = p0 - f(p0) * fp(p0) / (fp(p0)**2 - f(p0) * fpp(p0))
        # modified part iii
        p1 = p0- 4*f(p0)/fp(p0)

        p[it+1] = p1
        if (abs(p1-p0) < tol):
            pstar = p1
            info = 0
            return [p[:it+2],pstar,info,it]
        p0 = p1
    pstar = p1
    info = 1
    return [p,pstar,info,it]
            
driver()