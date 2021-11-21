import numpy as np


def e_vector(n, i):
    z = np.zeros((n, 1))
    z[i] = 1
    return z


def grad(f, x):
    n = len(x)
    m = len(f(x))
    dx = (abs(x) + 1) * np.finfo(float).eps ** (1.0 / 3.0)

    Df = np.zeros((n, m))
    for i in range(n):
        hi = e_vector(n, i) * dx
        df = 0.5 * (f(x + hi) - f(x - hi)) / hi[i]
        Df[i, :] = df.T

    return Df


def hessian(f, x):
    fx = f(x)
    n = len(x)
    dx = (abs(x) + 1) * np.finfo(float).eps ** 0.25
    H = np.zeros((n, n))
    for j in range(n):
        hj = e_vector(n, j) * dx
        for i in range(n):
            hi = e_vector(n, i) * dx
            H[i, j] = 0.25 * (f(x + hi + hj)
                              - f(x - hi + hj)
                              - f(x + hi - hj)
                              + f(x - hi - hj)
                              ) / (hi[i] * hj[j])
    return H
