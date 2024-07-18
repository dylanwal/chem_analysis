import numpy as np
from chem_analysis.utils.pad_edges import pad_edges_polynomial
from numpy.lib.stride_tricks import sliding_window_view

from chem_analysis.processing.processing_method import Smoothing
from chem_analysis.processing.weigths.weights import DataWeight


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
                  smoother: Smoothing = None,
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

    Returns
    -------
    mask:
        1 where the baseline is
        0 where peaks are
    """
    if not y.any():  # check if y is all zeros
        raise ValueError("y must have non-zero elements")
    if smoother is None:
        from chem_analysis.processing.smoothing.convolution import Gaussian
        smoother = Gaussian()

    window = max(window, 1)
    # compute noise level by breaking the data into sections and find section with min sigma
    stds = std_by_section(y, sections)
    min_sigma = np.percentile(stds, 5)  # min(stds)

    # smooth spectra with convolution
    _, smoothed_y = smoother.run(np.empty(0), y)

    # evaluate if point is outside min_sigma
    # (max(y_i) - min(y_i)) < n*sigma
    half_window = int(window / 2)
    padded_smoothed_y = pad_edges_polynomial(smoothed_y, degree=2, pad_amount=half_window)
    sliding_window = sliding_window_view(padded_smoothed_y, 2 * half_window + 1)
    mask = np.max(sliding_window, axis=1) - np.min(sliding_window, axis=1) < number_of_deviations * min_sigma

    if ignore_zeros:
        # added as other processing methods may set values to zero and should be ignored
        mask[y == 0] = False

    return mask



    min_sigma = np.percentile(y, 5 * number_of_deviations)
    # min_sigma = np.max(y) - abs(np.min(y))
    # for slice_ in slices:
    #     std_ = np.std(y[slice_])
    #     if std_ == 0:
    #         continue
    #     min_sigma = min(min_sigma, std_)


class MaxMinSigma(DataWeight):
    def __init__(self,
                 window: int = 3,
                 sections: int = 32,
                 number_of_deviations: int | float = 2,
                 smoother: Smoothing = None,
                 invert: bool = False,
                 ):
        """
        This algorithm assumes the y data contains at least one region with no signals.

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
