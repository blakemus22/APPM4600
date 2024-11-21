import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


def driver():
    N = 3 # degree of polynomial
    max_iter = 25 # maximum iterations of remez
    tol = 1e-6 # tolerance of remez algorithm

    f = lambda x: np.e**x 
    fp = lambda x: np.e**x

    pts = [-1,-.5,0,.5,1]

    (error,coef) = remez(pts,f, fp, max_iter,tol)
    print("coefficients for final model:", coef)
    print("error for final model", error)



def eval(x, coef):
    '''evaluates a function with coefficients coef at some value x'''

    res = 0
    for j in range(len(coef)):
        res+= coef[j] * x**j
    return res

def eval_deriv(x, coef):
    res = 0
    for j in range(1,len(coef)):
        res += coef[j]*j * x**(j-1)
    return res


def eval_error(f,x,coef):
   return f(x) - eval(x,coef)

def eval_error_deriv(fp, x, coef):
   return fp(x) - eval_deriv(x,coef)


def solve(nodes, fvals):
    '''takes the nodes and function values as input, outputs the coefficients for p(x) with
    oscillating error at those nodes'''
    A = np.zeros((len(fvals),len(fvals)))
    for i in range(len(fvals)):
        A[i][-1] = (-1)**i
        for j in range(len(fvals)-1):
            A[i][j] = nodes[i]**j

    Ainv = inv(A)

    x = inv(A) @ fvals
    return [x[:-1], x[-1]]

def bisection(f,func,a,b,tol,coefs,max_iter):
    '''input the endpoints and the function, find the root in that interval'''

    fa = f(func, a,coefs)
    fb = f(func, b, coefs)
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
      fd = f(func, d,coefs)
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

def remez(nodes, f, fp, max_iter, tol):
    '''calculates the remez approximation of the minimax approximation'''
    cur_nodes = nodes
    xvals = np.linspace(nodes[0],nodes[-1],500)

    for i in range(max_iter):
        fvals = np.array([f(cur_node) for cur_node in cur_nodes])
        # generate the coefficients of the approximation
        c, _ = solve(cur_nodes, fvals)
        print(c)

        # plot the error function
        funcvals = [eval_error(f,x,c) for x in xvals]
        plt.plot(xvals, funcvals, color = "indianred")
        plt.plot(cur_nodes, [eval_error(f,node,c) for node in cur_nodes], 'o', color = "slateblue")
        plt.vlines(cur_nodes,[0]*len(cur_nodes),[eval_error(f,node,c) for node in cur_nodes],linestyle = "-", color = "slateblue")
        plt.axhline(0, color = "black")
        plt.show()

        # generate the n+1 zeros of the function
        zeros = [bisection(eval_error, f,cur_nodes[j], cur_nodes[j+1], 1e-3, c, 25)[0] for j in range(len(cur_nodes)-1)]
        # now we need to find the maximum and minimum of the function between each of these zeros
        # generate the maximum error locations between these zeros
        new_nodes = [cur_nodes[0]] + [bisection(eval_error_deriv,fp, zeros[j],zeros[j+1], 1e-3, c, 25)[0] for j in range(len(zeros)-1)] + [cur_nodes[-1]]
        cur_nodes = new_nodes
        

        # check for stopping criteria
        error = [abs(eval_error(f,node,c)) for node in cur_nodes]
        if  max(error) - min(error) < tol:
            fvals = np.array([f(cur_node) for cur_node in cur_nodes])
            return (max(error),solve(cur_nodes,fvals)[0])


    # if we get here, the stopping criteria has been reached. not good
    fvals = np.array([f(cur_node) for cur_node in cur_nodes])   
    print("Hey just a heads up: we reached the maximum amount of iterations, you are probably in trouble!") 
    return (max(error),solve(cur_nodes,fvals)[0])


driver()