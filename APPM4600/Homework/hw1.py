#import the libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# problem 1

# get x values
x1vals = np.linspace(1.920,2.080,161)

# map onto y values
y_long = (x1vals ** 9) - 18*(x1vals ** 8) + 144*(x1vals ** 7) - 672*(x1vals ** 6) + 2016*(x1vals ** 5) - 4032*(x1vals ** 4) + 5376*(x1vals ** 3) - 4608*(x1vals ** 2) + 2304*(x1vals) - 512
y_short = (x1vals - 2)**9

# plot it all
plt.subplot(1,2,1)
plt.plot(x1vals, y_long)
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.title('Graph of (x-2)^9 - expanded version')

plt.subplot(1,2,2)
plt.plot(x1vals, y_short)
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.title('Graph of (x-2)^9 as a binomial to a power')

plt.show()

# problem 5

# create delta values
delta = 10**np.linspace(-16,0,17)

# map the delta values, using set values of x, to plot tie difference bewteen the functions

# part b
y_pi = abs((np.cos(np.pi + delta) - np.cos(np.pi)) + (2*np.sin(np.pi + delta/2)*np.sin(delta/2)))
y_big = abs((np.cos(10**6 + delta) - np.cos(10**6)) + (2*np.sin(10**6 + delta/2)*np.sin(delta/2)))

# part c
y_approx_pi = abs((np.cos(np.pi + delta) - np.cos(np.pi)) - (-1*delta*np.sin(np.pi)-delta**2/2*np.cos(np.pi + delta/2)))
y_approx = abs((np.cos(10**6 + delta) - np.cos(10**6)) - (-1*delta*np.sin(10**6)-delta**2/2*np.cos(10**6 + delta/2)))

plt.plot(np.log10(delta), y_approx_pi, label = "x = pi")
plt.plot(np.log10(delta), y_approx, label = "x = 10^6")

plt.xlabel('x axis: logarithmic scale : 10^x')
plt.ylabel('difference between expressions')
plt.legend()
plt.show()