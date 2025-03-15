
import numpy as np
from scipy.ndimage import gaussian_filter1d

from chem_analysis.utils.pad_edges import pad_edges_polynomial
from chem_analysis.processing.smoothing import Smoother
from numpy.lib.stride_tricks import sliding_window_view

from chem_analysis.processing.weigths.weights import Weights


def divide_array(array: np.ndarray, num_sections: int) -> list[slice]:
    section_length = len(array) // num_sections
    remainder = len(array) % num_sections

    slices = []
    start = 0
    for i in range(num_sections):
        if i < remainder:
            end = start + section_length + 1
        else:
            end = start + section_length
        slices.append(slice(start, end))
        start = end

    return slices


def std_by_section(array: np.ndarray, num_sections: int) -> np.ndarray:
    non_zero_array = array[array != 0]
    slices = divide_array(non_zero_array, num_sections)
    return np.array([np.std(non_zero_array[slice_]) for slice_ in slices])


def sectioned_std(y,
                  window: int = 3,
                  sections: int = 32,
                  number_of_deviations: int | float = 2,
                  smoother: Smoother | None = lambda y: gaussian_filter1d(y, 10),
                  ignore_zeros: bool = True,
                  ):
    """
    This algorithm assumes the y data contains at least one region with no signals.

    https://doi.org/10.1006/jmre.2000.2121

    Parameters
    ----------
    y:
        data
    window:
        number of points used to compute local variation
    sections:
        number of sections used to get minimum standard deviation || noise value
    number_of_deviations:
        the number of deviations from the noise value
        little effect, typically 2 to 4
    smoother:
        smoother used before doing min_max analysis
        highly recommended to use one. example: gaussian_filter1d(x, 10)
    ignore_zeros:
        won't apply mask if value is zero

    Returns
    -------
    mask:
        1 where the baseline is
        0 where peaks are
    """
    if not y.any():  # check if y is all zeros
        raise ValueError("y must have non-zero elements")

    window = max(window, 1)
    # compute noise level by breaking the data into sections and find section with min sigma
    stds = std_by_section(y, sections)
    min_sigma = np.percentile(stds, 5)  # min(stds)

    # smooth spectra
    smoothed_y = smoother(y) if smoother is not None else y

    # evaluate if point is outside min_sigma
    half_window = int(window / 2)
    padded_smoothed_y = pad_edges_polynomial(smoothed_y, degree=2, pad_amount=half_window)
    sliding_window = sliding_window_view(padded_smoothed_y, 2 * half_window + 1)
    mask = np.max(sliding_window, axis=1) - np.min(sliding_window, axis=1) < number_of_deviations * min_sigma

    if ignore_zeros:
        # added as other processing methods may set values to zero and should be ignored
        mask[y == 0] = False

    return mask


class SlidingWindowStd(Weights):
    def __init__(self,
                 window: int = 3,
                 sections: int = 32,
                 number_of_deviations: int | float = 2,
                 smoother: Smoother | None = lambda x: gaussian_filter1d(x, 10),
                 invert: bool = False,
                 ):
        """
        Finds the region with the smallest standard deviation.
        It uses the std to determine if a point within a window is outside that region.
        smoothing is usually applied to soften the analysis

        This algorithm assumes the y data contains at least one region with no signals which will be used to
        compute areas where variation exceeds the standard deviation.

        https://doi.org/10.1006/jmre.2000.2121

        Parameters
        ----------
        window:
            number of points used to compute local variation
        sections:
            number of sections used to get minimum standard deviation || noise value
        number_of_deviations:
            the number of deviations from the noise value
            little effect, typically 2 to 4
        smoother:
            smoother used before doing min_max analysis

        Returns
        -------
        mask:
            1 where the baseline is
            0 where peaks are
    """
        super().__init__(threshold=0.5, normalized=True, invert=invert)
        self.window = window
        self.sections = sections
        self.number_of_deviations = number_of_deviations
        self.smoother = smoother

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return sectioned_std(y, self.window, self.sections, self.number_of_deviations, self.smoother)
