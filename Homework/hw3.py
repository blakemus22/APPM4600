# import libraries
import numpy as np
import matplotlib.pyplot as plt

def driver():

# use routines    
    f = lambda x: x**3 + x -4
    a = 1
    b = 4

    f2 = lambda x: -np.sin(2*x) + 5*x/4 -3/4

    tol = 1e-10

    # bisection driver 

    [astar,ier, count] = bisection(f,a,b,tol)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('It took', count, 'iterations to get to the approximation')

    # fixed point driver

    x0 = 0
    Nmax = 100
    [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
    print('the approximate fixed point is:',xstar)
    print('f1(xstar):',f2(xstar))
    print('Error message reads:',ier)





# define routines

# used in questions 1, 2, and 3
# inputs were varied for each question, but code structure stayed the same
def bisection(f,a,b,tol):

    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,0]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,0]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,0]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier, count]
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
    return [astar, ier, count]

# used for question 5
def fixedpt(f,x0,tol,Nmax):


    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0)/abs(x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    
driver()

# problem 5

# define functions      
x_vals = np.linspace(-2,8, 1000)

f_5 = lambda x: x - 4 * np.sin(2*x) - 3 

conv_limit = lambda x: -np.sin(2*x) +5/4*x - 3/4

# plot functions for parts a and b
plt.plot(x_vals, f_5(x_vals), label = "function from part a")
plt.plot(x_vals, x_vals)
plt.plot(x_vals, conv_limit(x_vals), color = "green", label = "fixed point convergence")
plt.axhline(y = 0, color = "red", label = "y = 0")
plt.xlabel('x')
plt.ylabel('y')
plt.title("Finding a fixed point on the modified function")
plt.legend()
plt.show()

