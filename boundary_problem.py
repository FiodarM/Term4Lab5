__author__ = 'fiodar'

import numpy as np


def tridiag_solve(diags, f, conditions):
    a, b, c = diags
    alpha, beta = conditions
    x = np.zeros_like(f)
    k, l = np.zeros_like(f), np.zeros_like(f)

    k[0] = (alpha[2] - c[0] / c[1] * f[1]) / (a[0] - c[0] / c[1] * a[1])
    l[0] = (c[0] / c[1] - b[0]) / (a[0] - c[0] / c[1] * a[1])

    for i in xrange(1, len(f) - 1):
        k[i] = (f[i] - a[i] * k[i - 1]) / (a[i] * l[i - 1] + b[i])
        l[i] = - c[i] / (a[i] * l[i - 1] + b[i])

    x[-1] = (beta[2] - a[-1] * (k[-2] + l[-2] * k[-1]) - b[-1] * k[-1]) / \
            (a[-1] * l[-2] * l[-1] + b[-1] * l[-1] + c[-1])

    for i in xrange(1, len(f) - 1):
        x[-2 - i] = k[-2 - i] + l[-2 - i] * x[-1 - i]

    return x


def linear_2nd_order(coefs, bounds, conditions, n=50, f=lambda x: 0 * x):
    p, q = coefs
    x = np.linspace(bounds[0], bounds[1], n)
    h = x[1] - x[0]
    alpha, beta = conditions
    a = [alpha[0] / 2 - alpha[1] / h, 0.5 * beta[1] / h]
    b = [alpha[0] / 2 + alpha[1], -2 * beta[1] / h]
    c = [- 0.5 * alpha[1] / h, beta[0] + 0.5 * beta[1] / h]
    d = [alpha[2], beta[2]]

    a = np.insert(a, 1, 1 / h * (1 / h - 0.5 * p(x[1:-1])))
    b = np.insert(b, 1, q(x[1:-1]) - 2 / h ** 2)
    c = np.insert(c, 1, 1 / h * (1 / h + 0.5 * p(x[1:-1])))
    d = np.insert(d, 1, f(x[1:-1]))

    y = tridiag_solve((a, b, c), d, conditions)

    return y