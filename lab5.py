__author__ = 'fiodar'

from boundary_problem import *
import matplotlib.pyplot as plt

p = lambda x: 1. / x
q = lambda x: 1 - 0 / x ** 2
f = lambda x: -0. / x
a, b, n = 0., 3., 100
conditions = ([5., -4., -2.], [3., 1., 0.])
simple_conditions = ([1., 0, 0.], [1, 0, 1.])
y = linear_2nd_order((p, q), (a, b), simple_conditions, n, f)

plt.plot(np.linspace(a, b, n), y)
plt.show()