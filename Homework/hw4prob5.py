import math
import numpy as np
import matplotlib.pyplot as plt


def driver():


    x0 = 2

    tol = 1e-13

    f = lambda x : x**6 - x - 1
    fp = lambda x : 6*x**5 -1
    

# newtons method 
    [a,astar, ier, it] = newton(f,fp, x0, tol, 300)
    # print('the approximate root is', astar)
    # print('the error message reads:', ier)
    # print('f(astar)=', f(astar))
    # print(len(a), 'iterations to converge')
    err_lst_n = []
    for elt in range(len(a)-1):
        err_lst_n.append(abs(a[elt] - a[-1]))
    # print(err_lst_n)

    # plt.plot(err_lst_n[1:], err_lst_n[:-1])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylabel('x_k+1 - alpha')
    # plt.xlabel('x_k - alpha')
    # plt.title('Plotting consecutive errors on log-log axes: Newtons Method')
    # plt.show()

    n_slope = (math.log(err_lst_n[1]) - math.log(err_lst_n[-3]))/(math.log(err_lst_n[2]) - math.log(err_lst_n[-2]))
    print(n_slope)

# secant method 
    [elts, s, sier] = secant(f,x0,1,tol, 300)
    # print('the approximate root is', s)
    # print('the error message reads:', sier)
    err_lst_s = []
    for elt in range(len(elts)-1):
        err_lst_s.append(abs(elts[elt] - elts[-1]))
    # print(err_lst_s)

    plt.plot(err_lst_s[1:], err_lst_s[:-1], label = 'Secant convergence')
    plt.plot(err_lst_n[1:], err_lst_n[:-1], label = 'Newton convergence')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('x_k+1 - alpha')
    plt.xlabel('x_k - alpha')
    plt.title('Plotting errors on log-log axes: Newton/Secant Method')
    plt.legend()
    plt.show()

    s_slope = (math.log(err_lst_s[1]) - math.log(err_lst_s[-3]))/(math.log(err_lst_s[2]) - math.log(err_lst_s[-2]))
    print(s_slope)


# define routines

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

def secant(f, x0, x1, tol, Nmax):
    p = np.zeros(Nmax+2)
    p[0] = x0
    p[1] = x1
    if f(x0) == f(x1):
        return [p[:2], x1, 1]
    for j in range(Nmax):
        x2 = x1 - f(x1) * (x1-x0) / (f(x1) - f(x0))
        p[j+2] = x2
        if abs(x2 - x1) < tol:
            x2
            ier = 0
            return [p[:j+3], x2,ier]
        x0 = x1
        x1 = x2
        if f(x1) == f(x0):
            ier = 1
            return [p[:j+3], x1,ier]
    
    ier = 1
    return [p,x1,ier]
    
            
driver()