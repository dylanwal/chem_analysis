import numpy as np
from scipy import sparse

from chem_analysis.processing.processing_method import Baseline


def asymmetric_least_squared(
        y: np.ndarray,
        lambda_: float = 1e6,
        p: float = 1e-2,
        max_iter: int = 50,
) -> np.ndarray:
    """
   Asymmetric least squared (ALS) fitting.

    Parameters
    ----------
    y:
        y data
    lambda_:
        smoothing parameter
        larger values = smoother baselines
    p:
        penalizing weighting factor
        0 < p < 1
    max_iter:
        max number of fit iterations

    Returns
    -------
    baseline:

    """
    n = len(y)
    D = sparse.sparse.eye(n, format='csc')
    D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.
    D = D[1:] - D[:-1]
    D = D.T
    H = lambda_ * D.dot(D.T)
    w = np.ones(n)
    for i in range(max_iter):
        W = sparse.diags(w, 0, shape=(n, n))
        Z = W + H
        z = sparse.spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)

    return z


class AsymmetricLeastSquared(Baseline):
    def __init__(self,
                 lambda_: float = 1e6,
                 p: float = 1e-2,
                 max_iter: int = 50,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.lambda_ = lambda_
        self.p = p
        self.max_iter = max_iter

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        y_baseline, params = asymmetric_least_squared(
            y,
            self.lambda_,
            self.p,
            self.max_iter,
        )
        return y_baseline


def improved_asymmetric_least_squared(
        y: np.ndarray,
        lambda_: float = 1e6,
        p: float = 1e-2,
        max_iter: int = 50,
) -> np.ndarray:
    """
    Adaptive Iteratively Reweighted Penalized Least Squares (airPLS) baseline correction.


    Parameters
    ----------
    y:
        y data
    lambda_:
        Smoothness parameter.
        Higher values result in a smoother baseline.
    p:
       Order of the difference for the penalty (typically 2).
    max_iter:
        max number of fit iterations

    Returns
    -------
    baseline

    """
    # check inputs
    if max_iter < 2:
        raise ValueError("'max_iter' needs to be greater than 2")

    n = len(y)
    D = np.diff(np.eye(n), p, axis=0)  # Correct shape alignment
    H = lambda_ * D.T @ D  # Regularization term
    w = np.ones(n)

    for _ in range(max_iter):
        W = np.diag(w)
        Z = np.linalg.solve(W + H, w * y)  # Solve system
        d = y - Z
        w = np.exp(-1 * (d / (2 * (np.std(d) + 1e-8))) ** 2)  # Update weights with stability
        w /= np.max(w)  # Normalize weights

    return Z


class ImprovedAsymmetricLeastSquared(Baseline):
    def __init__(self,
                 lambda_: float = 1e6,
                 p: float = 2,
                 max_iter: int = 50,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.lambda_ = lambda_
        self.p = p
        self.max_iter = max_iter

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        y_baseline, params = improved_asymmetric_least_squared(
            y,
            self.lambda_,
            self.p,
            self.max_iter,
        )

        return y_baseline


# https://pubs.rsc.org/en/content/articlelanding/2010/an/b922045c
# https://diposit.ub.edu/dspace/bitstream/2445/188026/1/2014_IEEE_Adaptive_MarcoS_postprint.pdf  Adaptive Asymmetric Least Squares baseline estimation for analytical instruments     Sergio Oller-Moreno∗, Antonio Pardo‡
