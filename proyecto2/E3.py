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


def solve_problem(problem_provider, C, mu, params):
    n = C.shape[0]
    for param in params:
        f, h, g = problem_provider(C, mu, param)
        x, l, miu = solve(n, f, h, g)
        print('Solution for param: ', param)
        print('x:', x)
        print('l:', l)
        print('miu:', miu)
        print('f(x):', f(x))
        print('')


def main():
    # E3a
    n = 10
    C, mu = get_params(n)
    # E3b
    alfa_list = [np.min(mu), np.mean(mu)]
    print('Solving problem 1...')
    solve_problem(create_problem_1, C, mu, alfa_list)
    print('Done.\n')
    # E3c
    gamma_list = [1, 10, 100]
    print('Solving problem 3...')
    solve_problem(create_problem_3, C, mu, gamma_list)
    print('Done.\n')
    # E3d
    beta_list = [2, 1]
    print('Solving problem 2...')
    solve_problem(create_problem_2, C, mu, beta_list)
    print('Done.\n')


if __name__ == '__main__':
    main()
