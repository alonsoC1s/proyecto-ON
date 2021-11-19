import numpy as np


def get_params(
    n,
    diag_range=(0.01, 0.4),
    elem_range=(-0.2, 0.2),
    miu_range=(0.02, 0.1),
    random_state=123
):
    """Construye la matriz de covarianzas C y el
    vector de ganancia esperada de cada activo miu.
    Input: n
    Output: C, miu
    """
    rnp = np.random.RandomState(random_state)
    # Generate C
    elems = elem_range[0] + rnp.rand(n, n) * (elem_range[1] - elem_range[0])
    diag = diag_range[0] + rnp.rand(n, 1) * (diag_range[1] - diag_range[0])
    np.fill_diagonal(elems, diag)
    R = np.triu(elems)
    C = R.T @ R
    # Generate miu
    miu = miu_range[0] + rnp.rand(n, 1) * (miu_range[1] - miu_range[0])
    return C, miu


def find_min(f, x0):
    pass


def create_tolerance_generator():
    i = 1
    while True:
        yield 1.0 / i
        i += 1


def create_cost_generator():
    ck = 1
    while True:
        yield ck
        ck *= 2


def quadratic_penalty_method(
    Q,
    x0,
    c_generator=create_cost_generator(),
    tol_generator=create_tolerance_generator(),
    max_iter=1000
):
    xk = x0
    tolk = next(tol_generator)
    ck = next(c_generator)
    for k in range(max_iter):
        # a) Encontrar un mínimo x tal que minimice Qc iniciando en xk
        def Qc(x): return Q(x, ck)  # TODO: crear la Qc
        x = find_min(Qc, xk)
        # b) Ver si satisface las

        # c) si no lo hizo actualizar el costo y demás parámetros
        xk = x
        tolk = next(tol_generator)
        ck = next(c_generator)
    return xk


def create_Q(f, h, g):
    def Q(x, c):
        gg = g(x)
        gg[gg < 0] = 0
        return f(x) + 0.5*c*np.dot(h(x).T, h(x)) + 0.5*c*np.dot(gg.T, gg)
    return Q
    

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
    C, miu = get_params(n)
    f, h, g = create_problem_1(C, miu, 0.1)

    x0 = np.random.rand((n, 1))
    print('x:', x0)
    print('f(x):', f(x0))
    print('h(x):', h(x0))
    print('g(x):', g(x0))

    Q = create_Q(f, h, g)
    quadratic_penalty_method(Q, x0)