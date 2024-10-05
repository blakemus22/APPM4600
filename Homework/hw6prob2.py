import numpy as np
import math
from numpy.linalg import inv
from numpy.linalg import norm
def driver():
    x0 = np.array([-4,2,1])
    Nmax = 100
    tol = 1e-6

    [xstarN,ierN,itsN] = Newton(x0,tol,Nmax)
    print(xstarN)
    print('Newton: the error message reads: ',ierN)
    print('Newton: number of iterations is: ',itsN)

    [xstar,g,ier,its] = SteepestDescent(x0,tol,Nmax)
    print(xstar)
    print("Steepest Descent: the error messages reads: ", ier)
    print("Steepest Descent: number of iterations is: ", its)

    newx = SteepestDescent(x0,5e-2,Nmax)[0]
    firstit = SteepestDescent(x0,5e-2,Nmax)[3]
    print(firstit)
    [modx, modier, modits] = Newton(newx,tol,Nmax)
    print(modx)
    print("Hybrid Method: the error messages reads: ", modier)
    print("Hybrid Method: number of iterations is: ", modits+ firstit)




def evalF(x):
    F = np.zeros(3)
    F[0] = x[0] + np.cos(x[0]*x[1]*x[2])
    F[1] = (1-x[0])**.25 + x[1] + .05*x[2]**2 - .15*x[2] - 1
    F[2] = -x[0]**2 - .1*x[1]**2 + .01*x[1] + x[2] - 1
    return F
 
def evalJ(x):
    J = np.array([[1-x[1]*x[2]*np.sin(x[0]*x[1]*x[2]), -x[0]*x[2]*np.sin(x[0]*x[1]*x[2]), -x[0]*x[1]*np.sin(x[0]*x[1]*x[2])],
                  [.25/(1-x[0])**(.75), 1, .1*x[2] - .15],
                  [-2*x[0], -.2*x[1] + .01, 1]])
    return J

def eval_gradg(x):
    F = evalF(x)
    J = evalJ(x)
    
    gradg = np.transpose(J).dot(F)
    return gradg

def evalg(x):

    F = evalF(x)
    g = F[0]**2 + F[1]**2 + F[2]**2
    return g


def Newton(x0,tol,Nmax):

    for its in range(Nmax):
        J = evalJ(x0)
        Jinv = inv(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its+1]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]


###############################
### steepest descent code

def SteepestDescent(x,tol,Nmax):
    
    for its in range(Nmax):
        g1 = evalg(x)
        z = eval_gradg(x)
        z0 = norm(z)

        if z0 == 0:
            print("zero gradient")
        z = z/z0
        alpha1 = 0
        alpha3 = 1
        dif_vec = x - alpha3*z
        g3 = evalg(dif_vec)

        while g3>=g1:
            alpha3 = alpha3/2
            dif_vec = x - alpha3*z
            g3 = evalg(dif_vec)
            
        if alpha3<tol:
            print("no likely improvement")
            ier = 0
            return [x,g1,ier,its+1]
        
        alpha2 = alpha3/2
        dif_vec = x - alpha2*z
        g2 = evalg(dif_vec)

        h1 = (g2 - g1)/alpha2
        h2 = (g3-g2)/(alpha3-alpha2)
        h3 = (h2-h1)/alpha3

        alpha0 = 0.5*(alpha2 - h1/h3)
        dif_vec = x - alpha0*z
        g0 = evalg(dif_vec)

        if g0<=g3:
            alpha = alpha0
            gval = g0

        else:
            alpha = alpha3
            gval =g3

        x = x - alpha*z

        if abs(gval - g1)<tol:
            ier = 0
            return [x,gval,ier,its+1]

    print('max iterations exceeded')    
    ier = 1        
    return [x,g1,ier,Nmax]

driver()