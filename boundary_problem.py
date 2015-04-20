__author__ = 'fiodar'

import numpy as np


def tridiag_solve(diags, f, k0, l0):

    alpha, beta, gamma = diags
    x = np.zeros_like(f)
    k, l = np.zeros_like(f), np.zeros_like(f)
    k[0] = k0
    l[0] = l0
    for i in xrange(1, len(f) - 1):
        k[i] = (f[i] - alpha[i]*k[i-1])/(alpha[i]*l[i-1] + beta[i])
        l[i] = - gamma[i] / (alpha[i]*l[i-1] + beta[i])

    x[-1] = (f[-1] - alpha[-1]*k[-2]) / (alpha[-1]*l[-2] + beta[-1])
    for i in xrange(1, len(f) - 2):
        x[-2-i] = (f[-2-i] - alpha[-2-i]*k[-3-i] - gamma[-2-i]*x[-1-i]) /\
                  (alpha[-2-i]*l[-3-i] + beta[-2-i])

    x[0] = k[0] + l[0] * x[1]

    return x


def linear_2nd_order(coefs, f, bounds, conditions, n=50):
    p, q = coefs
    alpha, beta = conditions
    x = np.linspace(bounds[0], bounds[1], n)
    h = x[1] - x[0]

    a = [alpha[0] / 2 - alpha[1] / h, 0.5 * beta[1] / h]
    b = [alpha[0] / 2 + alpha[1], -2 * beta[1] / h]
    c = [- 0.5 * alpha[1] / h, beta[0] + 0.5 * beta[1] / h]
    d = [alpha[2], beta[2]]

    a = np.insert(a, 1, 1 / h * (1 / h - 0.5 * p(x[1:-1])))
    b = np.insert(b, 1, q(x[1:-1]) - 2 / h ** 2)
    c = np.insert(c, 1, 1 / h * (1 / h + 0.5 * p(x[1:-1])))
    d = np.insert(d, 1, f(x[1:-1]))


    k0 = 1 / (a[0] - c[0] / c[1] * a[1]) * f(x[1]) * (alpha[2] - c[0] / c[1] * f(x[1]))
    l0 = c[0] / c[1] - b[0]
    y = tridiag_solve((a, b, c), d, k0, l0)

    return y