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
plt.plot(x, p, color = "slateblue", label = "Minimax Approximation: $p(x) = 1.264 + 1.175x$")
plt.legend()
plt.title("Minimax approximation of $f(x) = e^x$ on the interval [-1,1]")

plt.subplot(1,2,2)
plt.plot(x,error, color = "slateblue")
plt.axhline(0, color = 'black')
plt.title("Error of the minimax approximation, $f(x) - p(x)$")
plt.vlines([-1,np.log(a1),1], [0,0,0], [np.e**(-1) - (a0-a1),a1 - (a1*np.log(a1)+a0), np.e - (a1+a0)], linestyle = "--", color = "indianred")
plt.plot([-1,np.log(a1),1],[np.e**(-1) - (a0-a1),a1 - (a1*np.log(a1)+a0), np.e - (a1+a0)], 'o', color = 'slateblue')
plt.show()



# x2 = np.linspace(.25,4,500)
# y = np.log(x2) + np.cos(x2) - x2/5
# plt.plot(x2,y, color = "indianred")
# plt.title("$f(x) = \ln(x) + \cos(x) - x/5$")
# plt.show()



error_vec = [.3599,.1029,.0351,.0208,.0113,.0059,.0032,.0017,.00095]
def compute_order(x):
    diff1 = np.abs(x[1::])

    diff2 = np.abs(x[0:-1])
    fit = np.polyfit(np.log(diff2.flatten()), np.log(diff1.flatten()), 1)
    print('the order of the equation is')
    print('log|(p_{n+1}-p}) = log(lambda) + alpha*log(|p_n-p|) where')
    print('lambda = ' + str(np.exp(fit[1])))
    print('alpha = ' + str(fit[0]))
    return [fit, diff1, diff2]

# compute_order(error_vec)