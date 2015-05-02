__author__ = 'fiodar'

from boundary_problem import *
import matplotlib.pyplot as plt

p = lambda x: 1. / x
q = lambda x: -3. * x
f = lambda x: -5. / x
a, b, n = 1., 3., 1000
conditions = ([5., -4., -2.], [3., 1., 0.])
y = ode_linear_2nd_order((p, q), (a, b), conditions, n, f)

plt.plot(np.linspace(a, b, n), y)
plt.show()