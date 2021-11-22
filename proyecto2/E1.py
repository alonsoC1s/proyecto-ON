import numpy as np


def get_params(
    n,
    diag_range=(0.01, 0.4),
    elem_range=(-0.2, 0.2),
    miu_range=(0.02, 0.1),
    random_state=123
):
    """Construye la matriz de covarianzas C y el
    vector de ganancia esperada de cada activo mu.
    Input: n
    Output: C, mu
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
