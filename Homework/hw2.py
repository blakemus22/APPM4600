import math
import matplotlib.pyplot as plt
import numpy as np
import random
# problem 2 
# part b

A = np.matrix([[.5,.5],[.5*(1+10**(-10)), .5*(1-10**(-10))]])
A_inv = np.matrix([[1-10**10, 10**10],[1+10**10, -10**10]])

print(np.linalg.norm(A, ord = 2), np.linalg.norm(A_inv, ord = 2))
print(np.linalg.eigvals(A))

# # part c
delta_b = np.matrix([10**-5, 2*10**-5])
delta_x = np.matmul(A_inv, np.transpose(delta_b))
print(delta_x)



# problem 3 

# part c

x = 9.999999995000000*10**-10
print(np.e**x - 1)

# part d
f = lambda x : x + x**2/2
print(f(x))


# problem 4

# part a
# initialize our vectors
# t = np.linspace(0, math.pi, 31)
# y = np.cos(t)
# s = 0
# for i in range(31):
#     s += t[i]*y[i]
# print("The sum is:", s)


# part b
R = 1.2
deltar = .1
f = 15
p = 0
theta = np.linspace(0, 2*math.pi, 400)
# x_theta = R * ((1 + deltar * np.sin(f*theta + p))*np.cos(theta))
# y_theta = R * ((1 + deltar * np.sin(f*theta + p))*np.sin(theta))


# plt.plot(x_theta, y_theta)
# plt.axis('equal')
# plt.show()

# part b
dr = .05

for i in range(1,11):
    R_2 = i
    f_2 = 2 + i
    p_2 = np.random.uniform(0,2)
    x_theta_loop = R_2 * ((1 + dr * np.sin(f_2*theta + p_2))*np.cos(theta))
    y_theta_loop = R_2 * ((1 + dr * np.sin(f_2*theta + p_2))*np.sin(theta))

    plt.plot(x_theta_loop, y_theta_loop)
plt.axis('equal')
plt.show()
