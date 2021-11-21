import numpy as np


def find_min(f, x0, tol):
    # TODO: implement
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
    f,
    h,
    g,
    x0,
    c_generator=create_cost_generator(),
    tol_generator=create_tolerance_generator(),
    max_iter=1000
):
    def Q(x, c):
        gg = g(x)
        gg[gg < 0] = 0
        return f(x) + 0.5*c*np.dot(h(x).T, h(x)) + 0.5*c*np.dot(gg.T, gg)

    xk = x0
    tolk = next(tol_generator)
    ck = next(c_generator)
    for k in range(max_iter):
        # a) Encontrar un mínimo x tal que minimice Qc iniciando en xk
        def Qc(x): return Q(x, ck)
        x = find_min(Qc, xk, tolk)
        # b) Ver si satisface las
        # TODO: Implement
        # c) si no lo hizo actualizar el costo y demás parámetros
        xk = x
        tolk = next(tol_generator)
        ck = next(c_generator)
    return xk
