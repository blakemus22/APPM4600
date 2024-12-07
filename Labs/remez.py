import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


def driver():
    N = 8 # degree of polynomial
    max_iter = 25 # maximum iterations of remez
    tol = 1e-9 # tolerance of remez algorithm

    # f = lambda x: np.log(x) + np.cos(x) - x/5
    # fp = lambda x: 1/x - np.sin(x) - 1/5

    f = lambda x: np.e**x 
    fp = lambda x: np.e**x

    a = -1
    b = 1
    pts_cheby = [a] + [b*np.cos((2*k-1)*np.pi/(2*N)) for k in range(1,N+1)][::-1] + [b]
    pts = [a + (b-a)/(N+1)*j for j in range(N+2)]

    (its,error,coef) = remez(pts,f, fp, max_iter,tol)
    print("coefficients for final model:", coef)
    print("error for final model", error)
    print("number of iterations:", its)

    xvals = np.linspace(pts[0],pts[-1],500)
    approx = [eval(x,coef) for x in xvals]
    fex = f(xvals)
    plt.plot(xvals, approx, color = "slateblue", label = "Approximation")
    plt.plot(xvals,fex,color = "indianred", label = "True function")
    plt.title("Final Approxmation")
    plt.legend()
    plt.show()



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
def special_bisection(f,funcp, func,a,b,tol, coefs,max_iter):
    '''this function is only to calculate the max value in the interval that touches the endpoints'''
    fpa = f(funcp,a,coefs)
    fpb = f(funcp,b,coefs)

    fa = f(func,a,coefs)
    fb = f(func,a,coefs)
    if (fpa*fpb >0):
        if abs(fa) > abs(fb):
            return [a,0]
        else:
            return [b,0]
    else:
        return bisection(f,func,a,b,tol,coefs,max_iter)


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

        # plot the error function
        funcvals = [eval_error(f,x,c) for x in xvals]
        plt.plot(xvals, funcvals, color = "indianred")
        plt.plot(cur_nodes, [eval_error(f,node,c) for node in cur_nodes], 'o', color = "slateblue")
        plt.vlines(cur_nodes,[0]*len(cur_nodes),[eval_error(f,node,c) for node in cur_nodes],linestyle = "--", color = "slateblue")
        plt.axhline(0, color = "black")
        plt.title(f'Remez Approximation at Iteration {i}')
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
            # final graph
            plt.plot(xvals, funcvals, color = "indianred")
            plt.plot(cur_nodes, [eval_error(f,node,c) for node in cur_nodes], 'o', color = "slateblue")
            plt.vlines(cur_nodes,[0]*len(cur_nodes),[eval_error(f,node,c) for node in cur_nodes],linestyle = "--", color = "slateblue")
            plt.axhline(0, color = "black")
            plt.title("Final Remez Approximation")
            plt.show()
        
            return (i+1,max(error),solve(cur_nodes,fvals)[0])


    # if we get here, the stopping criteria has not been reached. not good
    fvals = np.array([f(cur_node) for cur_node in cur_nodes])   
    print("Hey just a heads up: we reached the maximum amount of iterations, you are probably in trouble!") 
    return (max_iter, max(error),solve(cur_nodes,fvals)[0])


driver()