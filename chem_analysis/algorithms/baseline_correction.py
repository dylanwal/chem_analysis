from typing import Callable

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

BaselineCorrection = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]


def poly_baseline(x: np.ndarray, y: np.ndarray, deg: int = 0, mask: np.ndarray = None) \
        -> tuple[np.ndarray, np.ndarray]:

    if mask is not None:
        x_mask = x[np.where(mask)]
        y_mask = y[np.where(mask)]
    else:
        x_mask = x
        y_mask = y

    params = np.polyfit(x_mask, y_mask, deg)
    func_baseline = np.poly1d(params)
    y_baseline = func_baseline(x)
    y = y - y_baseline
    return x, y


def adaptive_polynomial_baseline(
        x: np.ndarray,
        y: np.ndarray,
        remove_amount: float = 0.4,
        deg: int = 0,
        num_iter: int = 5) \
        -> tuple[np.ndarray, np.ndarray]:
    """ Adaptive polynomial baseline correction

    The algorithm:
    1) calculate fit to the data
    2) removes some percent of data points furthest from the fit
    3) calculate a new fit from this 'masked data'
    4) starting with the raw data, removes some percent of data points (a larger percent then the first time) furthest
    from the fit
    5) repeat steps 3-4

    Parameters
    ----------
    x: np.ndarray[:]
        x data
    y: np.ndarray[:]
        y data
    remove_amount: float
        range: (0-1)
        amount of data to be removed by the end of all iterations
    deg: int
        range: [0, 10]
        degree of polynomial
    num_iter: int
        range: [2, +]
        number of iteration

    Returns
    -------

    """
    if num_iter < 1:
        raise ValueError("'num_iter' must be larger than 1.")

    x_mask = x
    y_mask = y
    number_of_points_to_remove_each_iteration = int(len(x) * remove_amount / num_iter)
    for i in range(num_iter):
        # perform fit
        params = np.polyfit(x_mask, y_mask, deg)
        func_baseline = np.poly1d(params)
        y_baseline = func_baseline(x)

        if i != num_iter:  # skip on last iteration
            # get values furthest from the baseline
            number_of_points_to_remove = number_of_points_to_remove_each_iteration * (i+1)
            index_of_points_to_remove = np.argsort(np.abs(y - y_baseline))[-number_of_points_to_remove:]
            y_mask = np.delete(y, index_of_points_to_remove)
            x_mask = np.delete(x, index_of_points_to_remove)

    y = y-y_baseline

    return x, y


def baseline_als(x: np.ndarray, y: np.ndarray, lam: float = 1_000, p: float = 0.01, niter=10):
    """Asymmetric Least Squares Smoothing """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return x, z
