__author__ = 'fiodar'

from boundary_problem import *
import matplotlib.pyplot as plt

p = lambda x: -4 * x ** 2
q = lambda x: -x ** 3
f = lambda x: - 3 * x
a, b, n = 1., 2., 500
conditions = ([5., 2., -1.], [2., 1., 3.])
y = ode_linear_2nd_order((p, q), (a, b), conditions, n, f)

plt.plot(np.linspace(a, b, n), y)
plt.xlim(1, 2)
plt.grid()
plt.show()