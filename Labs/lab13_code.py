import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time


def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     N = 5000
 
     ''' Right hand side'''
     b = np.random.rand(N,1)
     A = np.random.rand(N,N)
     built_in_start = time.time()
     x = scila.solve(A,b)
     built_in_end = time.time()
     print("Time for built in solver: " ,built_in_end - built_in_start)

     test = np.matmul(A,x)
     r = la.norm(test-b)
     
     # print(r)

     ## lu factorization
     lu_time_start = time.time()
     lu = scila.lu_factor(A)
     lu_time_end = time.time()
     lu_solve_start = time.time()
     x2 = scila.lu_solve(lu, b)
     lu_solve_end = time.time()
     print("Time for built in solver: " ,built_in_end - built_in_start)
     print("Time to compute LU factorization: ", lu_time_end - lu_time_start)
     print("Time to solve the LU factorization: ", lu_solve_end - lu_solve_start)

     ''' Create an ill-conditioned rectangular matrix '''
     N = 100
     M = 10
     A = create_rect(N,M)     
     b = np.random.rand(N,1)
     b[0] = .01

     # normal equations
     normal_solve_start = time.time()
     x1 = scila.solve(A.transpose() @ A, A.transpose() @ b)
     normal_solve_end = time.time()

     # QR method
     Q,R = la.qr(A)
     qr_start = time.time()
     x2 = scila.solve(R,Q.transpose() @ b)
     qr_end = time.time()
     print("error of normal equations: ", la.norm(x1 - b))
     print("error of qr equations: ", la.norm(x2 - b))
     print("Time to solve using normal equations: ", normal_solve_end - normal_solve_start)
     print("Time to solve using QR factorization: ", qr_end - qr_start)



     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,10,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
driver()       
