from numpy.linalg import norm
from .utils import grad
import numpy as np


def basic_interior_point(f, h, g, x0, s0, l0, miu0, gamma0, sigma, tau, tol = 1e-9, max_iter=1000):
    n = len(x0)
    p = len(l0)
    m = len(s0)

    def F(x, l, miu, s, gamma):
        gradf = grad(f, x)
        Ae = grad(h, x)
        Ai = grad(g, x)
        return np.concatenate((
            gradf + Ae.T@l + Ai.T@miu,
            h(x),
            g(x) + s,
            s * miu - gamma
        ))

    k = 0
    xk, sk, lk, miuk, gammak = x0, s0, l0, miu0, gamma0
    while norm(F(xk, lk, miuk, sk, 0)) > tol and k < max_iter:
        def Fk(x, l, miu, s): return F(x, s, l, miu, gammak)
        while norm(Fk(xk, lk, miuk, sk)) > gammak and k < max_iter:
            # Formamos el sistema (4.11)
            Jk = None # TODO: definirla
            # Resolvemos el sistema para obtener dirección de descenso
            d = np.linalg.solve(Jk, -Fk(xk, lk, miuk, sk))
            dx, dl, dm, ds = np.split(d, [n, p, m, m])
            # Calculamos alfas con (4.12)
            alfa_s = min(1, min(-tau * sk[ds < 0] / ds[ds < 0], default=1))
            alfa_miu = min(1, min(-tau*miuk[dm < 0] / dm[dm < 0], default=1))
            # Calcular nuevos valores
            k += 1
            xk += alfa_s * dx
            sk += alfa_s * ds
            lk += alfa_miu * dl
            miuk += alfa_miu * dm
            if k > max_iter:
                print('max iterations reached.')
    gammak *= sigma
    return xk, lk, miuk
        

def solve(n, f, h, g, random_state=42):
    rnp = np.random.RandomState(random_state)
    x0 = rnp.randn(n, 1)
    # To obtain lambda 0 we need Ae and Ai
    gradf0 = grad(f, x0)
    Ae = grad(h, x0)
    Ai = grad(g, x0)
    # Set gamma 0
    gamma0 = 1
    # Obtain s0 and miu0 that satisfies 4.10d
    m = Ae.shape[0]
    s0 = rnp.rand(m, 1)
    miu0 = gamma0 / s0
    # Obtain l0 with least-squares solution
    l0 = np.linalg.lstsq(Ae.T, -gradf0 -Ai.T@miu0)

    return basic_interior_point(
        f=f,
        h=h,
        g=g,
        x0=x0,
        s0=s0,
        l0=l0,
        miu0=miu0,
        gamma0=gamma0,
        sigma=rnp.rand(),
        tau=rnp.rand()
    )