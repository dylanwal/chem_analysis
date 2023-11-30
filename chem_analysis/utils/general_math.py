import logging

import numpy as np

MIN_FLOAT = np.finfo(float).eps


def pack_time_series(x: np.ndarray, time_: np.ndarray, z: np.array) -> np.ndarray:
    data = np.empty((len(time_) + 1, len(x) + 1), dtype=z.dtype)
    data[0, 0] = 0
    data[0, 1:] = x
    data[1:, 0] = time_
    data[1:, 1:] = z
    return data


def unpack_time_series(data: np.array) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[0, 1:]
    time_ = data[1:, 0]
    z = data[1:, 1:]
    return x, time_, z


def check_for_flip(x: np.ndarray, y: np.ndarray):
    """flip data if giving backwards; it should be low to high"""
    if x[0] > x[-1]:
        x = np.flip(x)
        y = np.flip(y)
    return x, y


def get_slice(x: np.ndarray, start=None, end=None, flip: bool = True) \
        -> slice:
    if start is None and end is None:
        return slice(0, len(x))
    elif start is None and end is not None:
        end_ = np.argmin(np.abs(x - end))
        return slice(0, end_)
    elif start is not None and end is None:
        start_ = np.argmin(np.abs(x - start))
        return slice(start_, len(x))

    if flip and start > end:
        start, end = end, start
    start_ = np.argmin(np.abs(x - start))
    end_ = np.argmin(np.abs(x - end))
    return slice(start_, end_)


def map_argmax_to_original(index: int | np.ndarray, mask) -> int | np.ndarray:
    """
    Map the index from masked array back to the index of the original array.

    Parameters:
    - argmax_index: Index obtained from the argmax of the masked array.
    - mask: Boolean mask indicating which elements were retained in the original array.

    Returns:
    - original_index: Index in the original array corresponding to the argmax of the masked array.
    """
    masked_indices = np.where(mask)[0]  # Get the indices of the retained elements
    original_index = masked_indices[index]
    return original_index


def normalize_by_max(y: np.ndarray) -> np.ndarray:
    return y / np.max(y)


def normalize_by_area(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x, y = check_for_flip(x, y)
    return y / np.trapz(x=x, y=y)


# pdf = probability distribution function
def get_mean_of_pdf(x: np.ndarray, y: np.ndarray = None, *, y_norm: np.ndarray = None) -> float:
    x, y = check_for_flip(x, y)
    if y_norm is None:
        y_norm = normalize_by_area(x, y)

    return np.trapz(x=x, y=x * y_norm)


def get_standard_deviation_of_pdf(x: np.ndarray, y: np.ndarray = None, *,
                                  y_norm: np.ndarray, mean: int | float = None
                                  ) -> float:
    x, y = check_for_flip(x, y)
    if y_norm is None:
        y_norm = normalize_by_area(x, y)
    if mean is None:
        mean = get_mean_of_pdf(x, y_norm=y_norm)

    return np.sqrt(np.trapz(x=x, y=y_norm * (x - mean) ** 2))


def get_skew_of_pdf(x: np.ndarray, y: np.ndarray = None, *,
                    y_norm: np.ndarray = None,
                    mean: int | float = None,
                    standard_deviation: int | float = None
                    ) -> float:
    x, y = check_for_flip(x, y)
    if y_norm is None:
        y_norm = normalize_by_area(x, y)
    if mean is None:
        mean = get_mean_of_pdf(x, y_norm=y_norm)
    if standard_deviation is None:
        standard_deviation = get_standard_deviation_of_pdf(x, y_norm=y_norm, mean=mean)

    return np.trapz(x=x, y=y_norm * (x - mean) ** 3) / standard_deviation ** 3


def get_kurtosis_of_pdf(x: np.ndarray, y: np.ndarray = None, *,
                        y_norm: np.ndarray = None,
                        mean: int | float = None,
                        standard_deviation: int | float = None
                        ) -> float:
    x, y = check_for_flip(x, y)
    if y_norm is None:
        y_norm = normalize_by_area(x, y)
    if mean is None:
        mean = get_mean_of_pdf(x, y_norm=y_norm)
    if standard_deviation is None:
        standard_deviation = get_standard_deviation_of_pdf(x, y_norm=y_norm, mean=mean)

    return (np.trapz(x=x, y=y_norm * (x - mean) ** 4) / standard_deviation ** 4) - 3


def get_full_width_at_height(x: np.ndarray, y: np.ndarray, height: float | int = 0.5) -> float:
    """ Calculates full width at a height. """
    lower, high = get_width_at(x, y, height)
    return abs(high - lower)


def get_asymmetry_factor(x: np.ndarray, y: np.ndarray, height: float | int = 0.1) -> float:
    """ Calculates asymmetry factor at height. """
    lower, high = get_width_at(x, y, height)
    middle = x[np.argmax(y)]

    return (high - middle) / (middle - lower)


def get_width_at(x: np.ndarray, y: np.ndarray, height: float | int = 0.5) -> tuple[float, float]:
    """ Determine full-width-x_max of a peaked set of points, x and y. """
    x, y = check_for_flip(x, y)
    height_half_max = np.max(y) * height
    index_max = np.argmax(y)
    if index_max == 0 or index_max == len(x):  # peak max is at end.
        logging.info("Finding fwhm is not possible with a peak max at an bound.")
        return 0, 0

    x_low = np.interp(height_half_max, y[:index_max], x[:index_max])
    x_high = np.interp(height_half_max, np.flip(y[index_max:]), np.flip(x[index_max:]))

    if x_low == x[0]:
        logging.info("fwhm or asym is having to linear interpolate on the lower end.")
        slice_ = max(3, int(index_max / 10))
        fit = np.polyfit(y[:slice_], x[:slice_], deg=1)
        p = np.poly1d(fit)
        x_low = p(height_half_max)

    if x_high == x[-1]:
        logging.info("fwhm or asym is having to linear interpolate on the lower end.")
        slice_ = max(3, int(index_max / 10))
        fit = np.polyfit(y[-slice_:], x[-slice_:], deg=1)
        p = np.poly1d(fit)
        x_high = p(height_half_max)

    return x_low, x_high
