import numpy as np
from scipy.ndimage import gaussian_filter1d, gaussian_filter, uniform_filter1d, uniform_filter
from scipy.signal.windows import gaussian

from chem_analysis.processing.processing_method import Smoothing


class Uniform(Smoothing):
    def __init__(self, size: float | int = 10, temporal_processing: int = 1):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.uniform_filter1d.html#scipy.ndimage.uniform_filter1d

        Parameters
        ----------
        size:
            Size of filter window
        """
        super().__init__(temporal_processing)
        self.size = size

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, uniform_filter1d(y, size=self.size)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, uniform_filter(z, size=self.size)

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


class Gaussian(Smoothing):
    def __init__(self, std: float | int = 10, temporal_processing: int = 1):
        """

        Parameters
        ----------
        std
            Standard deviation for Gaussian kernel.
        """
        super().__init__(temporal_processing)
        self.std = std

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, gaussian_filter1d(y, self.std)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, gaussian_filter(z, self.std)

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


def mollify(y: np.ndarray, kernel_size: int = 5, std: int | float = 2) -> np.ndarray:
    """
    Applies mollifier smoothing to a 1D signal.

    Parameters
    ----------
    y:
        The input signal.
    kernel_size:
        The size of the mollifier kernel.
        as kernal size increase it appraoches the result of 'Gaussian'
    std:
        The standard deviation of the Gaussian kernel.

    Returns
    -------

    """
    # Create a Gaussian kernel
    window = 2 * kernel_size + 1
    kernel = gaussian(window, std)
    kernel /= np.sum(kernel) 

    # Pad the signal to handle boundary effects
    padded_signal = np.pad(y, window // 2, mode='reflect')

    return np.convolve(padded_signal, kernel, mode='valid')


class Mollify(Smoothing):
    def __init__(self, kernel_size: int = 5, std: float | int = 10, temporal_processing: int = 1):
        """
        Applies mollifier smoothing to a 1D signal.
    
        Parameters
        ----------
        kernel_size:
            The size of the mollifier kernel.
            as kernal size increase it appraoches the result of 'Gaussian'
        std:
            The standard deviation of the Gaussian kernel.

        """
        super().__init__(temporal_processing)
        self.kernel_size = kernel_size
        self.std = std

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, mollify(y, self.kernel_size, self.std)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
