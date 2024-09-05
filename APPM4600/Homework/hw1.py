#import the libraries
import numpy as np
import matplotlib.pyplot as plt

#problem 1
x1vals = np.linspace(1.920,2.080,161)

y_long = (x1vals ** 9) - 18*(x1vals ** 8) + 144*(x1vals ** 7) - 672*(x1vals ** 6) + 2016*(x1vals ** 5) - 4032*(x1vals ** 4) + 5376*(x1vals ** 3) - 4608*(x1vals ** 2) + 2304*(x1vals) - 512
y_short = (x1vals - 2)**9

plt.plot(x1vals, y_long)
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.title('Graph of (x-2)^9 - expanded version')
plt.show()