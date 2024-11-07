import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(0,5,1000)

f1 = lambda x: x / (1+1/6*x**2 +7/360*x**4)
f2 = lambda x: (x - 7/60*x**3) / (1 + 1/20*x**2)
tay = lambda x: x - x**3/6 + x**5/120


est1 = f1(x)
est2 = f2(x)
taylor = tay(x)

fex = np.sin(x)


plt.semilogy(x, abs(est1 - fex), label = "part b")
plt.semilogy(x ,abs(est2 - fex), label = "part a and c")
plt.semilogy(x, abs(taylor - fex), label = "6th degree Taylor expansion")
plt.legend()
plt.title("Error using Pade approximations : sin(x)")
plt.show()