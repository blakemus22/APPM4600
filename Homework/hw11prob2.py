import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import quad


def driver():
    
    f = lambda x: x*np.cos(1/x)
    a = 1/100000000000
    b = 1


    SN = CompSimp(a,b,10,f)
    print("The approximate integral using Composite Simpsons rule is: ", SN)

   
    

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