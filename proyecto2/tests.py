import numpy as np
import unittest

import scipy.optimize

from E3 import create_problem_1, create_problem_2, create_problem_3
from E2 import solve
from E1 import get_params
from utils import e_vector, get_greatest_eigenvector, get_smallest_eigenvector


class ProblemasDeClase(unittest.TestCase):
    """
    Problemas de clase para probar la implementación
    """

    def test_problema_facil(self):
        """
            min x + y
            s.a     x^2 + y^2 = 2
        """
        def f(x): return np.array([[np.sum(x)]])
        def h(x): return np.array([[np.linalg.norm(x)**2 - 2]])
        def g(x): return np.zeros((1, 1))
        x, _, _ = solve(f, h, g, np.array([[-0.1], [-0.1]]))
        # x, _, _ = solve(f, h, g, np.array([[0.1], [0.1]])) # Se va al (1,1)
        print('\n', x)
        self.assertTrue(np.all(np.isclose(x, np.array([-1, -1]))))

    def test_problema_igualdades(self):
        n = 3
        G, _ = get_params(n, diag_range=(1.0, 2.0), elem_range=(1.0, 5.0), random_state=1)
        print('\nG:', G)
        def f(x): return -0.5 * x.T @ G @ x
        def h(x): return x.T @ x - 1
        def g(x): return -x
        sol = get_greatest_eigenvector(G)
        x, _, _ = solve(f, h, g, np.ones((n, 1)))
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol)))

    def test_problema_desigualdades(self):
        n = 3
        G, _ = get_params(n, diag_range=(1.0, 2.0), elem_range=(1.0, 2.0), random_state=2)
        print('\nG:', G)
        def f(x): return -0.5 * x.T @ G @ x
        def h(_): return np.zeros((1, 1))
        def g(x): return np.concatenate((x.T @ x - 1, -x * e_vector(n, 1)))
        sol = get_greatest_eigenvector(G)
        x, _, _ = solve(f, h, g, np.ones((n, 1)))
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol)))

    def test_problema_desigualdades_2(self):
        n = 3
        G, _ = get_params(n, diag_range=(1.0, 2.0), elem_range=(1.0, 2.0), random_state=123)
        def f(x): return 0.5 * x.T @ G @ x
        def h(_): return np.zeros((1, 1))
        def g(x): return np.concatenate((1 - x.T @ x, -x * e_vector(n, 1)))
        sol = get_smallest_eigenvector(G)
        x, _, _ = solve(f, h, g, np.array([[-0.76154248], [0.61654029], [-0.19982771]]))
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol)))

    def test_problem_1(self):
        # E3a
        n = 3
        C, mu = get_params(n)
        alfa = np.min(mu)
        f, h, g = create_problem_1(C, mu, alfa)
        x, _, _ = solve(f, h, g, np.ones((n, 1)))

        # get scipy solution
        cons = ({'type': 'eq', 'fun': h},
                {'type': 'ineq', 'fun': lambda x: -g(x)})
        res = scipy.optimize.minimize(f, np.ones((n, 1)), constraints=cons)
        sol = np.array([res.x]).T
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol, atol=0.01)))

    def test_problem_2(self):
        n = 3
        C, mu = get_params(n)
        beta = np.diag(C).min()
        f, h, g = create_problem_2(C, mu, beta)
        x, _, _ = solve(f, h, g, np.ones((n, 1)))

        # get scipy solution
        cons = ({'type': 'eq', 'fun': h},
                {'type': 'ineq', 'fun': lambda x: -(x.T @ C @ x - beta)},
                {'type': 'ineq', 'fun': lambda x: x})

        res = scipy.optimize.minimize(f, np.ones((n, 1)), constraints=cons)
        sol = np.array([res.x]).T
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol, atol=1e-5)))

    def test_problem_3(self):
        n = 3
        C, mu = get_params(n)
        gamma = 1
        f, h, g = create_problem_3(C, mu, gamma)
        x, _, _ = solve(f, h, g, np.ones((n, 1)))

        # get scipy solution
        cons = ({'type': 'eq', 'fun': h},
                {'type': 'ineq', 'fun': lambda x: -g(x)})
        res = scipy.optimize.minimize(f, np.ones((n, 1)), constraints=cons)
        sol = np.array([res.x]).T
        print('\nx:', x)
        print('\nf(x):', f(x))
        print('\nsol:', sol)
        print('\nf(sol):', f(sol))
        self.assertTrue(np.all(np.isclose(x, sol, atol=0.001)))


if __name__ == '__main__':
    unittest.main()
