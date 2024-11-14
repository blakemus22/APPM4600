import numpy as np
import matplotlib.pyplot as plt
import math

x = np.linspace(-1,1,500)
eee = np.e**x
guess = x+1.2

plt.plot(x,eee, color = "indianred", label = "$f(x) = e^x$")
plt.plot(x, guess, color = "slateblue", label = "Minimax guess, $f(x) = x + 1.2$")
plt.legend()
plt.title("Guess for minimax approximation")
plt.show()