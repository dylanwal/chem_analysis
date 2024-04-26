
import numpy as np
from scipy.ndimage import gaussian_filter

from chem_analysis.processing.processing_method import Smoothing


class ExponentialTime(Smoothing):
    def __init__(self, a: int | float = 0.8, non_temporal_processing: bool = False):
        """
        The exponential filter is a weighted combination of the previous estimate (output) with the newest input data,
        with the sum of the weights equal to 1 so that the output matches the input at steady state.

        Parameters
        ----------
        a:
            smoothing constant
            a is a constant between 0 and 1, normally between 0.8 and 0.99
        """
        super().__init__(non_temporal_processing)
        self.a = a
        self._other_a = 1-a

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError("Only valid for SignalArrays")

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


class GaussianTime(Smoothing):
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
        raise NotImplementedError("Only valid for SignalArrays")

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, gaussian_filter(z, self.sigma)

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
