import logging

import numpy as np

from chem_analysis.processing.processing_method import Baseline

logger = logging.getLogger(__name__)


def _LOESS(x: np.ndarray, y: np.ndarray, frac: float = 0.25) -> np.ndarray:
    """
    Locally weighted scatterplot smoothing (LOESS) implementation.

    Parameters
    ---------
    x :
    y:
    frac:
        The fraction of points used for smoothing
        0 < frac <= 1

    Returns
    -------
    baseline

    """
    y=np.copy(y)
    n = len(x)
    smoothed = np.zeros_like(y)

    # Distance matrix: Compute pairwise distances between points
    for i in range(n):
        # Compute the weight of the neighbors based on the distance
        distances = np.abs(x - x[i])
        bandwidth = np.percentile(distances, frac * 100)
        weights = np.exp(-0.5 * (distances / bandwidth)**2)  # Gaussian kernel

        # Weights for local regression
        X_local = np.vstack([np.ones_like(x), x - x[i]]).T  # Linear model: y = b0 + b1 * (x - x[i])
        W_local = np.diag(weights)

        # Compute the weighted least squares solution
        XTWX = X_local.T @ W_local @ X_local
        XTWY = X_local.T @ W_local @ y
        beta = np.linalg.solve(XTWX, XTWY)  # Solve for [b0, b1]

        # Compute the smoothed value for this point
        smoothed[i] = beta[0] + beta[1] * (x[i] - x[i])

    return smoothed


class LocallyWeightedScatterplotSmoothing(Baseline):
    def __init__(self,
                 frac: float = 0.25,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """


        Parameters
        ----------

    """
        super().__init__(temporal_processing, save_result)
        self.frac = frac

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return _LOESS(x, y, self.frac)


LOESS = LocallyWeightedScatterplotSmoothing


def _RLOESS(x: np.ndarray, y: np.ndarray, frac: float = 0.25, max_iter: int = 10) -> np.ndarray:
    """
    Robust Locally Weighted Scatterplot Smoothing (LOESS) with iterative reweighting

    Parameters
    ---------
    x :
    y:
    frac:
        The fraction of points used for smoothing
        0 < frac <= 1
    max_iter:
        max iteration

    Returns
    -------
    baseline

    """
    n = len(x)
    r = int(np.ceil(frac * n))

    # Initial weights
    weights = np.ones(n)

    for _ in range(max_iter):
        y_est = np.zeros(n)

        for i in range(n):
            # Distances and local weights
            distances = np.abs(x - x[i])
            sorted_indices = np.argsort(distances)

            # Take r closest points
            neighborhood_indices = sorted_indices[:r]

            # Calculate tri-cube weights
            max_distance = np.max(distances[neighborhood_indices])
            local_weights = (1 - (distances[neighborhood_indices] / max_distance)**3)**3 if max_distance > 0 else np.ones(len(neighborhood_indices))

            # Apply robustness weights
            combined_weights = weights[neighborhood_indices] * local_weights

            # Weighted least squares
            A = np.vstack([np.ones(r), x[neighborhood_indices]]).T
            b = y[neighborhood_indices]

            #Solve the equation
            theta = np.linalg.solve(A.T @ np.diag(combined_weights) @ A, A.T @ np.diag(combined_weights) @ b)
            y_est[i] = theta[0] + theta[1] * x[i]

        # Calculate residuals and update robustness weights
        residuals = y - y_est

        # Median absolute deviation of residuals
        mad = np.median(np.abs(residuals - np.median(residuals)))

        # Robustness weights using bisquare function
        s = residuals / (6 * mad)
        weights = (1 - s**2)**2 * (np.abs(s) < 1)

    return y_est


class RobustLocallyWeightedScatterplotSmoothing(Baseline):
    def __init__(self,
                 frac: float = 0.25,
                 max_iter: int = 10,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """


        Parameters
        ----------

    """
        super().__init__(temporal_processing, save_result)
        self.frac = frac
        self.max_iter = max_iter

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return _RLOESS(x, y, self.frac, self.max_iter)


RLOESS = RobustLocallyWeightedScatterplotSmoothing
