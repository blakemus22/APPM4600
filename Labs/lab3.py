# import libraries
import numpy as np

def driver_bisection():

# use routines 
    # problem 1
    # f = lambda x: x**2*(x-1)
    # a = -1
    # b = 2

    # problem 2 
    f = lambda x: np.sin(x)
    a = .5
    b = 3*np.pi/4


    tol = 1e-5

    [astar,ier] = bisection(f,a,b,tol)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))




def driver_fixedpt():

# test functions 
     f1 = lambda x: x - (x**5 - 7) /12


     Nmax = 100
     tol = 1e-10

# test f1 '''
     x0 = 1
     [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
    



# define routines
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b)
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
      fd = f(d)
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
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier]
      

# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
# verify that these 4 functions have 7^(1/5) as a fixed point
# f1 = lambda x: x*(1+(7-x**5)/x**2)**3
# f2 = lambda x: x - (x**5 - 7) / x**2
# f3 = lambda x: x - (x**5 - 7) / (5*x**4) 
# f4 = lambda x: x - (x**5 - 7) /12
# print(f1(7**(1/5)), f2(7**(1/5)), f3(7**(1/5)), f4(7**(1/5)))
driver_fixedpt()
            
