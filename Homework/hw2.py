import math
import matplotlib.pyplot as plt
import numpy as np



# problem 4

# part a
# initialize our vectors
t = np.linspace(0, math.pi, 31)
y = np.cos(t)
# this part is very confusing. I'm not sure what the notation means

# part b
R = 1.2
deltar = .1
f = 15
p = 0
theta = np.linspace(0, 2*math.pi, 400)
x_theta = R * (1 + deltar * np.sin(f*theta + p)*np.cos(theta))
y_theta = R * (1 + deltar * np.sin(f*theta + p)*np.sin(theta))



#plt.plot(x_theta, y_theta)
#plt.axis('equal')
#plt.show()

# part b
dr = .05

for i in range(1,11):
    R_2 = i
    f_2 = 2 + i
    p_2 = np.random.uniform(0,2)
    x_theta_loop = R_2 * (1 + dr * np.sin(f_2*theta + p_2)*np.cos(theta))
    y_theta_loop = R_2 * (1 + dr * np.sin(f_2*theta + p_2)*np.sin(theta))

    plt.plot(x_theta_loop, y_theta_loop)
plt.axis('equal')
plt.show()
