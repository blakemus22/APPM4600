import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import quad


def driver():
    
    f = lambda x: 1 / (1+x**2)
    a = -5
    b = 5
    
    # exact integral
    #I_ex = i'm not sure I would have to do some math, but I could figure it out


    # nvals = [n for n in range(20,201,10)]
    # trap_val = [CompTrap(a,b,n,f) for n in range(20,201,10)]
    # simp_val = [CompSimp(a,b,n,f) for n in range(20,201,10)]

    # plt.plot(nvals, trap_val, color = "pink", label = "Composite Trapezoidal Approximation")
    # plt.plot(nvals, simp_val, color = "blue", label = "Composite Simpson Approximation")
    # plt.legend()
    # plt.title("Approximations of the integral $\int_{-5}^5 1/(1+x^2)dx$")
    # plt.xlabel("Number of nodes")
    # plt.show()


    # part c: using scipy quadrature

    TN = CompTrap(a,b,1292,f)

    SN = CompSimp(a,b,96,f)

    smalltol = quad(f,a,b,full_output = 1, epsabs = 10e-6)
    biggertol = quad(f,a,b,full_output = 1, epsabs = 10e-4)

    print("Number of Function evaluations to obtain TN :", 1292)
    print("Number of Function evaluations to obtain SN :", 96)
    print("Number of Function evaluations for scipy with tolerance 10e-4: ", biggertol[2]['neval'])
    print("Number of Function evaluations for scipy with tolerance 10e-6: ", smalltol[2]['neval'])   
   

  
def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    
    I_trap = h*f(xnode[0])*1/2
    
    for j in range(1,n):
         I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    
    return I_trap     

def CompSimp(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_simp = f(xnode[0])

    nhalf = n/2
    for j in range(1,int(nhalf)+1):
         # even part 
         I_simp = I_simp+2*f(xnode[2*j])
         # odd part
         I_simp = I_simp +4*f(xnode[2*j-1])
    I_simp= I_simp + f(xnode[n])
    
    I_simp = h/3*I_simp
    
    return I_simp     


    
    
driver() 