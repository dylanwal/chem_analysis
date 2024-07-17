
import numpy as np
from scipy.signal import savgol_filter

from chem_analysis.processing.processing_method import Smoothing


class SavitzkyGolay(Smoothing):
    """

    """
    def __init__(self, window_length: int = 10, order: int = 3, temporal_processing: int = 1):
        """
        The Savitzky Golay filter is a particular type of low-pass filter, well adapted for data smoothing.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.

        Parameters
        ----------
        window_length:
            The length of the filter window (i.e., the number of coefficients)
            window_length must be less than or equal to the size of y
        order:
            The order of the polynomial used to fit the samples.
            order must be less than window_length.
        """
        super().__init__(temporal_processing)
        if order > window_length:
            raise ValueError(f"'SavitzkyGolay.order'({order}) must be less than window_length ({window_length}).")
        self.window_length = window_length
        self.order = order

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.window_length > len(y):
            raise ValueError(f"'SavitzkyGolay.window_length'({self.window_length}) must be less than or "
                             f"equal to the size of y ({len(y)}).")
        return x, savgol_filter(y, self.window_length, self.order)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.window_length > z.shape[0]:
            raise ValueError(f"'SavitzkyGolay.window_length'({self.window_length}) must be less than or "
                             f"equal to the first dimension of z ({z.shape[0]}).")
        return x, y, savgol_filter(z, self.window_length, self.order, axis=0)

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    # def run_2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    #     # TODO: double check if this truely does 2D
    #     if self.window_length > len(z.shape[0]) or self.window_length > len(z.shape[1]):
    #         raise ValueError(f"'SavitzkyGolay.window_length'({self.window_length}) must be less than or "
    #                          f"equal to both dimensions of z ({len(z.shape)}).")
    #     return x, y, savgol_filter(z, self.window_length, self.order)








