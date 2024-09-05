x = [1,2,3]

import numpy as np 
y = np.array([1,2,3])
import matplotlib.pyplot as plt
X = np.linspace(0,2 * np.pi, 100)
Ya = np.sin(X)
Yb = np.cos(X)


x = np.linspace(0,10,100)
y = np.arange(100)

print("the first three entries of x are", x[:3])


w = 10**(-np.linspace(1,10,10))
print(w)

x_new = np.linspace(1,10,10)
s = 3*w

plt.plot(X, Ya)
plt.plot(X, Yb)
plt.plot(x_new,s)
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.show()