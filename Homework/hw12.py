import math
from numpy.linalg import inv
from numpy.linalg import norm
import numpy as np

def driver():
    # part a

    # for n in range(4,21,4):
    #     mat = create_h(n)
    #     eigen, error ,its = power_method(mat,1e-5, n,50)
    #     print("The dominant eigenvalue in the", n,"x",n, "matrix:", eigen)
    #     print("Number of iterations:", its)
    #     print("Error message:", error)

    # part b

    mat = create_h(16)
    eigen, _, _ = power_method(mat,1e-5, 16,50)
    lambdaI = np.diag(np.full(16,eigen))
    newmat = mat - lambdaI
    print(eigen)
    eigen2, error2 ,its2 = power_method(newmat,1e-2, 16,100)
    print(eigen2)
    print("The smallest eigenvalue in the 16 x 16 matrix:", eigen2 + eigen)
    print("Number of iterations:", its2)
    print("Error message:", error2)






def power_method(mat, tol,n, max_it):
    init = np.ones(n)
    x0 = mat @ init
    for i in range(max_it):
        x1 = mat @ x0
        diff = x0/x0[-1] - x1/x1[-1]
        if norm(diff) < tol:
            # we found the eigenvector
            x1 = x1 / find_mag(x1)

            e_val = x1.dot(mat @ x1)
            return [e_val, 0, i]
        else:
            x0 = x1


    # if we get here, we never found an eigenvector
    return [0, 1, max_it] 

def find_mag(vec):
    mag = 0
    for el in range(len(vec)):
        mag += vec[el]**2

    mag = math.sqrt(mag)
    return mag


def create_h(n):
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i][j] = 1 / (i+j+1)

    return A

driver()