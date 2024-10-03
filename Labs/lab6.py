import math
import numpy as np


# pre lab right here !! 
h = .01 * 1 / 2 ** (np.arange(0,10))

forward_diff = np.cos(np.pi/2 + h) / h
# print(forward_diff)

center_diff = (np.cos(np.pi/2 + h) - np.cos(np.pi/2 - h)) / (2*h)

err = (center_diff[1:] + 1 )/ (center_diff[:-1] + 1)
# print(center_diff)
# print(err)

# both of these converge very quickly(same rate). I think that in this case, it is linear convergence. 
# the error reduction between each step betweens to have a constant ratio.
# I calculated the ratio between consecutive errors and we get about .25 for each step. Therefore, it is linear.


# real lab

import numpy as np
import math
from numpy.linalg import inv
from numpy.linalg import norm
def driver():
    x0 = np.array([1,0])
    Nmax = 100
    tol = 1e-10
    [xstar,ier,its] = Newton(x0,tol,Nmax)
    print(xstar)
    print('Newton: the error message reads:',ier)
    print('Netwon: number of iterations is:',its)

    # [xstar, ier, its] = SlackerNewton(x0, tol, 1e-2, 100)
    # print(xstar)
    # print('Lazy Newton: the error message reads:',ier)
    # print('Lazy Newton: number of iterations is:',its)


def evalF(x):
    F = np.zeros(2)
    F[0] = 4*x[0]**2 + x[1]**2 - 4
    F[1] = x[0] + x[1] - np.sin(x[0] - x[1])
    return F

def evalJ(x, h):
    # this is pretty confusing code, but I am trying to use forward difference here
    J = np.array([[(4*(x[0]+h)**2 - 4*(x[0] +h)**2)/h, (x[1]+h)**2 - (x[1])**2], 
                  [(((h - np.sin(x[0] +h - x[1])) + np.sin(x[0] - x[1])))/h,
                    (((h - np.sin(x[0] +h - x[1])) + np.sin(x[0] - x[1])))/h]])
    return J


def Newton(x0,tol,Nmax):

    for its in range(Nmax):
        J = evalJ(x0, 10**(-3)* abs(x0))
        Jinv = inv(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def SlackerNewton(x0,tol,tolupdate, Nmax):

    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
        elif norm(x1-x0) < tolupdate:
            x0 = x1
            J = evalJ(x0)
            Jinv = inv(J)

        else:
            x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

driver()