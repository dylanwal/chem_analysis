import numpy as np
from scipy.ndimage import gaussian_filter1d, gaussian_filter, uniform_filter1d, uniform_filter

from chem_analysis.processing.processing_method import Smoothing


class Uniform(Smoothing):
    def __init__(self, size: float | int = 10, non_temporal_processing: bool = False):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.uniform_filter1d.html#scipy.ndimage.uniform_filter1d

        Parameters
        ----------
        size
            Size of filter window
        """
        super().__init__(non_temporal_processing)
        self.size = size

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, uniform_filter1d(y, size=self.size)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, uniform_filter(z, size=self.size)

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


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
        return x, gaussian_filter1d(y, self.sigma)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, gaussian_filter(z, self.sigma)

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
