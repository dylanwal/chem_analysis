import abc
from typing import Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.masks import DataMask


class BaselineCorrection(ProcessingMethod, abc.ABC):
    def __init__(self, mask: DataMask | Iterable[DataMask] = None):
        self.mask = mask

    def _apply_mask(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.mask is None:
            return x, y
        if isinstance(self.mask, DataMask):
            return self.mask.get_data(x, y)

        for mask in self.mask:
            x, y = mask.get_data(x, y)

        return x, y

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class Polynomial(BaselineCorrection):
    def __init__(self, degree: int = 1, mask: DataMask | Iterable[DataMask] = None):
        super().__init__(mask)
        self.degree = degree
        self.mask = mask
        self._y_baseline = None

    @property
    def y_baseline(self) -> np.ndarray:
        return self._y_baseline

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x, y = self._apply_mask(x, y)

        params = np.polyfit(x, y, self.degree)
        func_baseline = np.poly1d(params)
        y_baseline = func_baseline(x)
        y = y - y_baseline

        self._y_baseline = y_baseline
        return x, y

# Whittaker smoother (filter=21 Hz, smoothfactor = 32768)
# Bernstein polynomial (order = 3)
# Ablative (points 5, passes 10)
# Splines
# multipoint

from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix, eye, diags


def asls(data, lam=1e6, p=1e-2, diff_order=2, max_iter=50, tol=1e-3, weights=None):
    """
    Fits the baseline using assymetric least squared (AsLS) fitting.

    Parameters
    ----------
    data : array-like, shape (N,)
        The y-values of the measured data, with N data points.
    lam : float, optional
        The smoothing parameter. Larger values will create smoother baselines.
        Default is 1e6.
    p : float, optional
        The penalizing weighting factor. Must be between 0 and 1. Residuals
        above the data will be given p weight, and residuals below the data
        will be given p-1 weight. Default is 1e-2.
    diff_order : {2, 1, 3}, optional
        The order of the differential matrix. Default is 2 (second order
        differential matrix).
    max_iter : int, optional
        The max number of fit iterations. Default is 50.
    tol : float, optional
        The exit criteria. Default is 1e-3.
    weights : array-like, shape (N,), optional
        The weighting array. If None (default), then the initial weights
        will be an array with size equal to N and all values set to 1.

    Returns
    -------
    z : numpy.ndarray, shape (N,)
        The calculated baseline.
    params : dict
        A dictionary with the following items:

        * 'weights': numpy.ndarray, shape (N,)
            The weight array used for fitting the data.

    Raises
    ------
    ValueError
        Raised if p is not between 0 and 1.

    Notes
    -----
    Algorithm initially developed in [1]_ and [2]_, and code was adapted from
    https://stackoverflow.com/questions/29156532/python-baseline-correction-library.

    References
    ----------
    .. [1] Eilers, P. A Perfect Smoother. Analytical Chemistry, 2003, 75(14),
           3631-3636.
    .. [2] Eilers, P., et al. Baseline correction with asymmetric least squares
           smoothing. Leiden University Medical Centre Report, 2005, 1(1).

    """
    if p < 0 or p > 1:
        raise ValueError('p must be between 0 and 1')
    y, D, W, w = _setup_whittaker(data, lam, diff_order, weights)
    for _ in range(max_iter):
        z = spsolve(W + D, w * y)
        mask = (y > z)
        w_new = p * mask + (1 - p) * (~mask)
        if relative_difference(w, w_new) < tol:
            break
        w = w_new
        W.setdiag(w)

    residual = y - z
    return z, {'roughness': z.T * D * z, 'fidelity': residual.T * W * residual, 'weights': w}

def relative_difference(a, b):
    # Implement this function or replace it with your existing implementation
    return np.sum(np.abs(a - b)) / np.sum(np.abs(a + b))


def _setup_whittaker(data, lam, diff_order=2, weights=None):
    """
    Sets the starting parameters for doing penalized least squares.

    Parameters
    ----------
    data : array-like, shape (M,)
        The y-values of the measured data, with M data points.
    lam : float
        The smoothing parameter, lambda. Typical values are between 10 and
        1e8, but it strongly depends on the penalized least square method
        and the differential order.
    diff_order : {2, 1, 3, 4, 5}, optional
        The integer differential order; either 1, 2, 3, 4, or 5. Default is 2.
    weights : array-like, shape (M,), optional
        The weighting array. If None (default), then will be an array with
        shape (M,) and all values set to 1.

    Returns
    -------
    y : numpy.ndarray, shape (N,)
        The y-values of the measured data, converted to a numpy array.
    scipy.sparse.dia.dia_matrix
        The product of lam * D.T * D, where D is the sparse diagonal matrix of
        the differential, and D.T is the transpose of D.
    scipy.sparse.dia.dia_matrix
        The sparse weight matrix with the weighting array as the diagonal values.
    weight_array : numpy.ndarray, shape (N,), optional
        The weighting array.

    """
    y = np.asarray(data)
    diff_matrix = difference_matrix(y.shape[0], diff_order)
    if weights is None:
        weight_array = np.ones(y.shape[0])
    else:
        weight_array = np.asarray(weights).copy()
    return y, lam * diff_matrix.T * diff_matrix, diags(weight_array), weight_array


def difference_matrix(data_size, diff_order=2):
    """
    Creates an n-order differential matrix.

    Parameters
    ----------
    data_size : int
        The number of data points.
    diff_order : {2, 1, 3, 4, 5}, optional
        The integer differential order; either 1, 2, 3, 4, or 5. Default is 2.

    Returns
    -------
    scipy.sparse.dia.dia_matrix
        The sparse diagonal matrix of the differential.

    Raises
    ------
    ValueError
        Raised if diff_order is not 1, 2, 3, 4, or 5.

    Notes
    -----
    Most baseline algorithms use 2nd order differential matrices when
    doing penalized least squared fitting.

    It would be possible to support any differential order by doing
    np.diff(np.eye(data_size), diff_order), but the resulting matrix could
    cause issues if data_size is large. Therefore, it's better to only
    provide sparse arrays for the most commonly used differential orders.

    The resulting matrices are transposes of the result of
    np.diff(np.eye(data_size), diff_order). Not sure why there is a discrepancy,
    but this implementation allows using the differential matrices are they
    are written in various publications, ie. D.T * D rather than having to
    do D * D.T like most code, such as those adapted from stack overflow:
    https://stackoverflow.com/questions/29156532/python-baseline-correction-library.

    """
    if diff_order not in (1, 2, 3, 4, 5):
        raise ValueError('The differential order must be 1, 2, 3, 4, or 5')
    diagonals = {
        1: [-1, 1],
        2: [1, -2, 1],
        3: [-1, 3, -3, 1],
        4: [1, -4, 6, -4, 1],
        5: [-1, 5, -10, 10, -5, 1]
    }[diff_order]

    return diags(diagonals, list(range(diff_order + 1)), shape=(data_size - diff_order, data_size))