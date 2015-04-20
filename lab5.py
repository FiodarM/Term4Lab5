__author__ = 'fiodar'

from boundary_problem import *
import matplotlib.pyplot as plt

p = lambda x: 1. / x
q = lambda x: -3 * x
f = lambda x: -5. / x
a, b = 1, 3
conditions = ([5., -4., -2.], [3., 1., 0.])
y = linear_2nd_order((p, q), f, (a, b), conditions)

plt.plot(np.linspace(a, b), y)
plt.show()