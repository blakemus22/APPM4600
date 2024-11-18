import numpy as np
import matplotlib.pyplot as plt
import math

x = np.linspace(-1,1,500)
eee = np.e**x
guess = x+1.2

# plt.plot(x,eee, color = "indianred", label = "$f(x) = e^x$")
# plt.plot(x, guess, color = "slateblue", label = "Minimax guess, $p(x) = x + 1.2$")
# plt.legend()
# plt.title("Guess for minimax approximation")
# plt.show()

a1 = (np.e - 1/np.e)/2
a0 = (np.e - a1*np.log(a1))/2
p = a1 * x + a0

error = eee - p
plt.subplot(1,2,1)
plt.plot(x,eee, color = "indianred", label = "$f(x) = e^x$")
plt.plot(x, p, color = "slateblue", label = "p(x):  Minimax approximation")
plt.legend()
plt.title("Calculated minimax approximation")

plt.subplot(1,2,2)
plt.plot(x,abs(error), color = "slateblue")
plt.title("Absolute error of the approximation")
plt.show()
