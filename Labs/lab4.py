# import libraries
import numpy as np
import math
    
def driver():

# test functions 
     f1 = lambda x: (10/(x+4))**.5


     Nmax = 100
     tol = 1e-10

# test f1 '''
     x0 = 1.5
     out = fixedpt(f1,x0,tol,Nmax)
#      print('the approximate fixed point is:',out[-1])
#      print(out)
#     # print('f1(xstar):',f1(out[-1]))
#     # print('Error message reads:',ier)

#      compute_order(out, 1.3652300134140976)
     print(aitkens(out, 10, 10))
    



# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    X = np.zeros((Nmax,1))
    count = 0
    while (count <Nmax):
       X[count] = x0
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          X[count] = xstar
          ier = 0
          return X[:count+1]
       x0 = x1

    xstar = x1
    X[count + 1] = xstar
    ier = 1
    return X
    

def compute_order(x, xstar):
    diff1 = np.abs(x[1::] - xstar)

    diff2 = np.abs(x[0:-1]-xstar)

    fit = np.polyfit(np.log(diff2.flatten()), np.log(diff1.flatten()), 1)
    print('the order of the equation is')
    print('log|(p_{n+1}-p}) = log(lambda) + alpha*log(|p_n-p|) where')
    print('lambda = ' + str(np.exp(fit[1])))
    print('alpha = ' + str(fit[0]))
    return [fit, diff1, diff2]


def aitkens(x, tol, max):
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    return x0 - (x1-x0)**2/(x2 -2*x1 + x0)

driver()
