import numpy as np
from E1 import get_params
from E2 import solve


def create_problem_1(C, miu, alfa):
    def f(x): return x.T @ C @ x
    def h(x): return np.ones(x.shape).T @ x - 1

    def g(x):
        g1 = alfa - miu.T @ x
        g2 = -x
        return np.concatenate((g1, g2))
    return f, h, g


def create_problem_2(C, miu, beta):
    def f(x): return miu.T @ x
    def h(x): return np.ones(x.shape).T @ x - 1

    def g(x):
        g1 = x.T @ C @ x - beta
        g2 = -x
        return np.concatenate((g1, g2))
    return f, h, g


def create_problem_3(C, miu, gamma):
    def f(x): return -miu.T @ x + gamma * x.T @ C @ x
    def h(x): return np.ones(x.shape).T @ x - 1
    def g(x): return -x
    return f, h, g


if __name__ == '__main__':
    n = 3
    C, mu = get_params(n)
    alfa = np.min(mu)
    f, h, g = create_problem_1(C, mu, alfa)
    x, l, miu = solve(n, f, h, g)
    print('x:', x)
    print('l:', l)
    print('miu:', miu)
