import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


def driver():
    N = 3 # degree of polynomial
    max_iter = 25 # maximum iterations of remez
    tol = 1e-3 # tolerance of remez algorithm

    f = lambda x: np.e**(2*x) * np.sin(x)
    fp = lambda x: np.e**(2*x) * np.cos(x) + 2* np.e**(2*x) * np.sin(x)

    






def eval(x, coef):
    '''evaluates a function with coefficients coef at some value x'''

    res = 0
    for j in range(len(coef)):
        res+= coef[j] * x**j
    return res

def eval_deriv(x, coef):
    res = 0
    for j in range(1,len(coef)):
        res += coef[j]*j * x**j-1


def eval_error(f,x,coef):
   return f(x) - eval(x,coef)

def eval_error_deriv(fp, x, coef):
   return fp(x) - eval_deriv(x,coef)


def solve(nodes, fvals):
    '''takes the nodes and function values as input, outputs the coefficients for p(x) with
    oscillating error at those nodes'''
    A = np.zeros(len(fvals),len(fvals))
    for i in range(len(fvals)):
        A[i][-1] = (-1)**j
        for j in range(len(fvals)-1):
            A[i][j] = nodes[i]**j

    Ainv = inv(A)

    x = inv(A) @ fvals
    return [x[:-1], x[-1]]

def bisection(f,a,b,tol,coefs,max_iter):
    '''input the endpoints and the function, find the root in that interval'''

    fa = f(a,coefs)
    fb = f(b, coefs)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d,coefs)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count +=1
      if count > max_iter:
         return [d,1]
      
    astar = d
    ier = 0
    return [astar, ier]

def remez(nodes, f, max_iter, tol):
    '''calculates the remez approximation of the minimax approximation'''
    cur_nodes = nodes

    for i in range(max_iter):
        fvals = np.array([f(cur_node) for cur_node in cur_nodes])
        # generate the coefficients of the approximation
        c, _ = solve(cur_nodes, fvals)

        # generate the n+1 zeros of the function
        zeros = [bisection(eval_error, cur_nodes[j], cur_nodes[j+1], 1e-3, c, 25) for j in range(len(cur_nodes)-1)]
        # now we need to find the maximum and minimum of the function between each of these zeros

        # generate the maximum error locations between these zeros
        new_nodes = [cur_nodes[0]] + [bisection(eval_error_deriv,zeros[j],zeros[j+1], 1e-3, c, 25) for j in range(len(zeros)-1)] + [cur_nodes[-1]]
        cur_nodes = new_nodes
        # check for stopping criteria
        if max(cur_nodes) - min(cur_nodes) < tol:
            fvals = np.array([f(cur_node) for cur_node in cur_nodes])
            return solve(cur_nodes,fvals)[1]


    # if we get here, the stopping criteria has been reached. not good
    fvals = np.array([f(cur_node) for cur_node in cur_nodes])   
    print("Hey just a heads up: we reached the maximum amount of iterations, you are probably in trouble!") 
    return solve(cur_nodes,fvals)