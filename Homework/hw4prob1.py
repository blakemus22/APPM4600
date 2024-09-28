import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import special


#plotting the temperature function at x feet deep after 60 days:
alpha = .136e-6
xvals = np.linspace(0,5,1000)
yvals = 35* scipy.special.erf(xvals / (2*math.sqrt(alpha*5184000))) - 15

plt.plot(xvals, yvals)
plt.axhline(y = 0, color = 'green')
plt.xlabel("Depth of water line")
plt.ylabel("Temperature after 60 days")
plt.title("Temperature of a pipeline based on depth over time")

plt.show()
def driver():

    alpha = .138e-6
    t = 5184000
# use routines    
    # problem 1
    f = lambda x: 35* scipy.special.erf(x / (2*math.sqrt(alpha*t))) - 15
    fp = lambda x: 70 / (2*math.sqrt(alpha*t)*math.sqrt(np.pi)) * np.e**(-1*(x/(2*math.sqrt(alpha*t)))**2)


    a = 0
    b = 5
    x0 = 5

    tol = 1e-13
    

    [astar,ier] = bisection(f,a,b,tol)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))

    [a, astar, ier, it] = newton(f,fp, x0, tol, 200)
    print('the approximate root is', astar)
    print('the error message reads:', ier)
    print('f(astar)=', f(astar))
    print(it, 'iterations to converge')


# define routines
def bisection(f,a,b,tol): 

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
      
    astar = d
    ier = 0
    return [astar, ier]

def newton(f,fp, p0,tol,Nmax):

    p = np.zeros(Nmax+1)
    p[0] = p0
    for it in range(Nmax):
        p1 = p0- f(p0)/fp(p0)

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