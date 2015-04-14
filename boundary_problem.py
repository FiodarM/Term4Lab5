__author__ = 'fiodar'

import numpy as np

def tridiag_solve(diags, f):

    alpha, beta, gamma = diags
    x = np.zeros_like(f)
    k, l = np.zeros_like(f), np.zeros_like(f)
    k[0] = f[0] / beta[0]
    l[0] = - gamma[0] / beta[0]
    for i in range(1, len(f) - 1):
        k[i] = (f[i] - alpha[i]*k[i-1])/(alpha[i]*l[i-1] + beta[i])
        l[i] = - gamma[i] / (alpha[i]*l[i-1] + beta[i])

    x[-1] = (f[-1] - alpha[-1]*k[-2]) / (alpha[-1]*l[-2] + beta[-1])
    for i, item in enumerate(x[-2:0:-1]):
        x[-2-i] = (f[-2-i] - alpha[-2-i]*k[-3-i] - gamma[-2-i]*x[-1-i]) /\
                  (alpha[-2-i]*l[-3-i] + beta[-2-i])

    x[0] = k[0] + l[0] * x[1]

    return x

