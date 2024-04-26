import numpy as np
from scipy.ndimage import gaussian_filter

from chem_analysis.processing.processing_method import Smoothing


class Gaussian(Smoothing):
    def __init__(self, sigma: float | int = 10, non_temporal_processing: bool = False):
        """

        Parameters
        ----------
        sigma
            Standard deviation for Gaussian kernel.
        """
        super().__init__(non_temporal_processing)
        self.sigma = sigma

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, gaussian_filter(y, self.sigma)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, gaussian_filter(z, self.sigma)

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
