import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math

def driver():


    f = lambda x: 1/(1+x**2)
    fp = lambda x: -2*x/(1+x**2)**2

    N = 20
    ''' interval'''
    a = -5
    b = 5
   
    ''' create equispaced interpolation nodes'''
    xint = [5*np.cos((2*j+1)*np.pi/(2*(N+1))) for j in range(N+1)]
    xint = xint[::-1]
    
    ''' create interpolation data'''
    yint = np.zeros(N+1)
    ypint = np.zeros(N+1)
    for jj in range(N+1):
        yint[jj] = f(xint[jj])
        ypint[jj] = fp(xint[jj])
    
    ''' create points for evaluating the Lagrange interpolating polynomial'''
    Neval = 1000
    xeval = np.linspace(a,b,Neval+1)
    yevalL = np.zeros(Neval+1)
    yevalH = np.zeros(Neval+1)
    for kk in range(Neval+1):
      yevalL[kk] = eval_lagrange(xeval[kk],xint,yint,N)
      yevalH[kk] = eval_hermite(xeval[kk],xint,yint,ypint,N)

    '''natural cubic spline'''
    (M,C,D) = create_natural_spline(yint,xint,N)
    yevalNC = eval_cubic_spline(xeval,Neval,xint,N,M,C,D)

    '''clamped cubic spline'''
    (M2,C2,D2) = create_clamped_spline(ypint, yint, xint, N)
    yevalCC = eval_cubic_spline(xeval, Neval, xint, N, M2, C2, D2)

    ''' create vector with exact values'''
    fex = np.zeros(Neval+1)
    for kk in range(Neval+1):
        fex[kk] = f(xeval[kk])

    
    plt.subplot(1,2,1)
    plt.plot(xeval,fex,'ro-', label = 'True Function', ms = 1)
    plt.plot(xeval,yevalL,'bs--',label='Lagrange', ms = .8) 
    plt.plot(xeval,yevalH,'c--',label='Hermite', ms = .8)
    plt.plot(xeval,yevalNC,'g--', label = 'Natural Cubic Spline', ms = .8)
    plt.plot(xeval,yevalCC,'y--', label = 'Clamped Cubic Spline', ms = .8)
    plt.legend()
    plt.title("Interpolation using Chebyshev Nodes")
         
    errL = abs(yevalL-fex)
    errH = abs(yevalH-fex)
    errNC = abs(yevalNC - fex)
    errCC = abs(yevalCC - fex)
    plt.subplot(1,2,2)
    plt.semilogy(xeval,errL,'bs--',label='Lagrange', ms = .8)
    plt.semilogy(xeval,errH,'c--',label='Hermite', ms = .8)
    plt.semilogy(xeval,errNC,'g--', label = 'Natural Cubic Spline', ms = .8)
    plt.semilogy(xeval,errCC, 'y--', label = 'Clamped Cubic Spline', ms = .8)
    plt.legend()
    plt.title("Error Using Chebyshev Nodes: N = 20")
    plt.show()            


def eval_hermite(xeval,xint,yint,ypint,N):

    ''' Evaluate all Lagrange polynomials'''

    lj = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    ''' Construct the l_j'(x_j)'''
    lpj = np.zeros(N+1)
#    lpj2 = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
#              lpj2[count] = lpj2[count]*(xint[count] - xint[jj])
              lpj[count] = lpj[count]+ 1./(xint[count] - xint[jj])
              

    yeval = 0.
    
    for jj in range(N+1):
       Qj = (1.-2.*(xeval-xint[jj])*lpj[jj])*lj[jj]**2
       Rj = (xeval-xint[jj])*lj[jj]**2

       yeval = yeval + yint[jj]*Qj+ypint[jj]*Rj
       
    return(yeval)
       


def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)
  
def create_natural_spline(yint,xint,N):

#    create the right  hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)
    for i in range(1,N):
       hi = xint[i] - xint[i-1]
       hip = xint[i+1] - xint[i]

       b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
       h[i-1] = hi
       h[i] = hip
    h[N] = xint[N] - xint[N-1]

#  create the matrix A so you can solve for the M values
    A = np.zeros((N+1,N+1))
    A[0][0] = 1
    A[N][N] = 1
    for r in range(1,N):
        A[r][r-1] = h[r-1]/6
        A[r][r] = (h[r-1] + h[r])/3
        A[r][r+1] = h[r]/6
            
#  Invert A    
    Ainv = la.inv(A)

# solver for M    
    M  = la.inv(A) @ b
    
#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = yint[j]/h[j] - h[j]*M[j]/6
       # find the C coefficients
       D[j] = yint[j+1]/h[j] - h[j]*M[j+1]/6
       # find the D coefficients
    return(M,C,D)


def create_clamped_spline(ypint, yint,xint,N):
#    create the right  hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N)
   #  vector values

    for i in range(1,N):
       hi = xint[i] - xint[i-1]
       hip = xint[i+1] - xint[i]

       b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
       h[i-1] = hi
       h[i] = hip

    b[0] = -1*ypint[0] + (yint[1] - yint[0])/(h[0])
    b[-1] = -1*ypint[-1] + (yint[-1] - yint[-2])/h[N-1]

#  create the matrix A so you can solve for the M values
    A = np.zeros((N+1,N+1))
    # first and last rows of equation
    A[0][0] = h[0]/3
    A[0][1] = h[0]/6
    A[N][N] = h[-1]/3
    A[N][N-1] = h[-1]/6
    # create the rest of the matrix
    for r in range(1,N):
        A[r][r-1] = h[r-1]/6
        A[r][r] = (h[r-1] + h[r])/3
        A[r][r+1] = h[r]/6
            

#  Invert A    
    Ainv = la.inv(A)

# solver for M    
    M  = la.inv(A) @ b
    
#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = yint[j]/h[j] - h[j]*M[j]/6
       # find the C coefficients
       D[j] = yint[j+1]/h[j] - h[j]*M[j+1]/6
       # find the D coefficients
    
    return(M,C,D)
       
def eval_local_spline(xeval,xi,xip,Mi,Mip,C,D):
# Evaluates the local spline as defined in class
# xip = x_{i+1}; xi = x_i
# Mip = M_{i+1}; Mi = M_i

    hi = xip-xi
   
    yeval = (Mi / (hi*6))*(xip - xeval)**3 + (Mip/(hi*6)) *(xeval - xi)**3 + C*(xip-xeval) + D*(xeval-xi)
    return yeval 
    
    
def  eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    
    yeval = np.zeros(Neval+1)
    
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        
#   find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]

# evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
#   copy into yeval
        yeval[ind] = yloc
    ind1 = np.where((xeval < xint[0]))
    xloc1 = xeval[ind1]
    yloc1 = eval_local_spline(xloc1,xint[0], xint[1],M[1],M[2],C[1],D[1])
    yeval[ind1] = yloc1

    indlast = np.where((xeval > xint[-1]))
    xloclast = xeval[indlast]
    yloclast = eval_local_spline(xloclast,xint[-2], xint[-1],M[-2],M[-1],C[-1],D[-1])
    yeval[indlast] = yloclast


    return(yeval)

           
driver()   

     