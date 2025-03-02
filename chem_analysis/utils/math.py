import logging
import copy
from typing import Sequence

import numpy as np

logger = logging.getLogger(__name__)

MIN_FLOAT = np.float64.min  # np.finfo(float).eps
DTYPE_UINT = (np.uint8, np.uint16, np.uint32, np.uint64)
DTYPE_INT = (np.int8, np.int16, np.int32, np.int64)
DTYPE_FLOAT = (np.float16, np.float32, np.float64)


def norm2(old: np.ndarray, new: np.ndarray) -> float:
    if np.all(old == 0):
        return 0.0
    return np.linalg.norm(new - old) / np.linalg.norm(old)


def pack_time_series(x: np.ndarray, time_: np.ndarray, z: np.array) -> np.ndarray:
    if x.shape[0] != z.shape[1]:
        raise ValueError(f"'x.shape[0]' must equal 'z.shape[1]'\n\tx shape:{x.shape}\n\tz shape:{z.shape}")
    if time_.shape[0] != z.shape[0]:
        raise ValueError(f"'time_.shape[0]' must equal 'z.shape[0]'\n\ttime shape:{time_.shape}\n\tz shape:{z.shape}")

    data = np.empty((len(time_) + 1, len(x) + 1), dtype=z.dtype)
    data[0, 0] = 0
    data[0, 1:] = x
    data[1:, 0] = time_
    data[1:, 1:] = z
    return data


def unpack_signal2D(data: np.array) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def quick_check_for_sorted_array(x: np.ndarray, min_check: int = 5000) -> bool:
    """
    ** Not a strict check ** but it is quick to compute for any size array
    Parameters
    ----------
    x
    min_check

    Returns
    -------

    """
    if len(x) < 2:
        return True
    if x[0] > x[-1]:
        return False

    if len(x) < min_check:
        return np.all(np.all(x[:-1] <= x[1:]))

    i = np.random.randint(1, len(x)-1, min_check-2)
    i.sort()
    return np.all(x[i[:-1]] <= x[i[1:]])


def get_slice(
        x: np.ndarray,
        start: int | float | None = None,
        end: int | float | None = None,
        *,
        checks: bool = True,
        strict_bounds: bool = True,
        start_bound: bool | None = None,
        end_bound: bool | None = None,
) -> slice:
    """
    gets slice from the nearest values

    Parameters
    ----------
    x
        sorted list small -> big
    start:
        value to start slice
    end:
        value to end slice
    checks:
        checks for issues with inputs
    strict_bounds:
        True: raises ValueError if bounds can't be satisfied
        False: puts None in slice
    start_bound:
        None: closest
        False: lowest
        True: highest
    end_bound:
        None: closest
        False: lowest
        True: highest
    Returns
    -------

    """
    if start is None and end is None:
        return slice(None, None)

    if checks:
        if start is not None and end is not None and start > end:
            raise ValueError("'start' value is larger than 'end'. \nFix: Flip bounds.")
        if not quick_check_for_sorted_array(x):
            raise ValueError("Array is not sorted. \nFix: sort 'x'")
        if start is not None and start < np.min(x):
            raise ValueError(f"'start' is outside the range of 'x'. \n\tstart:{start}\n\tmin(x): {np.min(x)}")
        if end is not None and end > np.max(x):
            raise ValueError(f"'end' is outside the range of 'x'. \n\tend:{end}\n\tmax(x): {np.max(x)}")

    if start is None:
        start_ = None
    elif start_bound is None:
        start_ = np.argmin(np.abs(x - start))
    else:
        if start_bound:
            if end is None:
                mask = (x >= start)
            else:
                mask = (end >= x) & (x >= start)
        else:
            mask = (x <= start)
        if np.all(mask == 0):
            if strict_bounds:
                raise ValueError("slice can't find value for start.")
            else:
                start_ = None
        else:
            start_ = np.argmin(np.abs(x[mask] - start)) + np.argmax(mask)

    if end is None:
        end_ = None
    elif end_bound is None:
        end_ = np.argmin(np.abs(x - end))
    else:
        if end_bound:
            mask = (x >= end)
        else:
            mask = (x <= end) & (x >= start)
        if np.all(mask == 0):
            if strict_bounds:
                raise ValueError("slice can't find value for start.")
            else:
                end_ = None
        else:
            end_ = np.argmin(np.abs(x[mask] - end)) + np.argmax(mask) + 1

    return slice(start_, end_)


def get_slice_by_nearest_y(y: np.ndarray, value: int | float) -> slice:
    """"""
    max_index = np.argmax(y)
    if max_index == 0:
        bound = np.argmin(np.abs(y - value))
        return slice(None, bound)
    if max_index == len(y) - 1:
        bound = np.argmin(np.abs(y - value))
        return slice(bound, None)

    lb = np.argmin(np.abs(y[:max_index] - value))
    ub = np.argmin(np.abs(y[max_index:] - value)) + max_index
    return slice(lb, ub)


def map_argmax_to_original(index: int | np.ndarray, mask) -> int | np.ndarray:
    """
    Map the index from a masked array back to the index of the original array.

    Parameters
    ----------
    index : int
        Index obtained from the argmax of the masked array.
    mask : array_like, bool
        Boolean mask indicating which elements were retained in the original array.

    Returns
    -------
    int or np.ndarray
        Index in the original array corresponding to the argmax of the masked array.
        Returns None if argmax_index is out of bounds or the corresponding mask value is False.
    """
    masked_indices = np.where(mask)[0]  # Get the indices of the retained elements
    original_index = masked_indices[index]
    return original_index


def get_index_of_values_in_common(array_a: np.ndarray, array_b: np.ndarray) -> list[int | None]:
    """ get index where values in array_b is in array_a """
    indices = []
    for value in array_b:
        result = np.where(array_a == value)[0]
        indices.append(result[0] if len(result) > 0 else None)

    return indices


def normalize_by_max(y: np.ndarray) -> np.ndarray:
    return y/np.max(y)


def normalize_by_max_with_x_range(y: np.ndarray, x: np.ndarray, x_range: Sequence[int | float] = None) -> np.ndarray:
    if x_range[0] > x_range[1]:
        x_range = copy.copy(x_range)
        x_range = x_range[1], x_range[0]

    slice_ = get_slice(x, *x_range)
    return y / np.max(y[slice_])


def normalize_by_area(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x, y = check_for_flip(x, y)
    return y / np.trapz(x=x, y=y)


def y_normalized_by_area_with_x_range(x: np.ndarray, y: np.ndarray, x_range: Sequence[int | float] = None) \
        -> np.ndarray:
    if x_range[0] > x_range[1]:
        x_range = copy.copy(x_range)
        x_range = x_range[1], x_range[0]

    slice_ = get_slice(x, *x_range)
    return y / np.trapz(x=x[slice_], y=y[slice_])


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
    """
    skew: float
    skew
    symmetric: -0.5 to 0.5; moderate skew: -1 to -0.5 or 0.5 to 1; high skew: <-1 or >1;
    positive tailing to higher numbers; negative tailing to smaller numbers

    Returns
    -------

    """
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
    """
        kurtosis: float
            kurtosis (Fisher) (Warning: highly sensitive to peak bounds)
            negative: flatter peak; positive: sharp peak

        Returns
        -------

    """
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
    """
    asymmetry_factor: float
    asymmetry factor; distance from the center line of the peak to the back slope divided by the distance from the
    center line of the peak to the front slope;
    >1 tailing to larger values; <1 tailing to smaller number

    Returns
    -------

    """
    lower, high = get_width_at(x, y, height)
    middle = x[np.argmax(y)]

    return float((high - middle) / (middle - lower))


def get_width_at(x: np.ndarray, y: np.ndarray, height: float | int = 0.5) -> tuple[float, float]:
    """ Determine full-width-x_max of a peaked set of points, x and y. """
    x, y = check_for_flip(x, y)
    height_half_max = np.max(y) * height
    index_max = np.argmax(y)
    if index_max == 0 or index_max == len(x):  # peak max is at end.
        logger.info("Finding fwhm is not possible with a peak max at an bound.")
        return 0, 0

    x_low = np.interp(height_half_max, y[:index_max], x[:index_max])
    x_high = np.interp(height_half_max, np.flip(y[index_max:]), np.flip(x[index_max:]))

    if x_low == x[0]:
        logger.info("fwhm or asym is having to linear interpolate on the lower end.")
        slice_ = max(3, int(index_max / 10))
        fit = np.polyfit(y[:slice_], x[:slice_], deg=1)
        p = np.poly1d(fit)
        x_low = p(height_half_max)

    if x_high == x[-1]:
        logger.info("fwhm or asym is having to linear interpolate on the lower end.")
        slice_ = max(3, int(index_max / 10))
        fit = np.polyfit(y[-slice_:], x[-slice_:], deg=1)
        p = np.poly1d(fit)
        x_high = p(height_half_max)

    return x_low, x_high


def min_uint_dtype(max_value: int | float) -> np.dtype:
    # Find the smallest dtype that can represent all values in the array
    for dtype in DTYPE_UINT:
        if max_value <= np.iinfo(dtype).max:
            return dtype

    raise ValueError(f"Number is larger than '{DTYPE_UINT[-1]}' and can't be represented.")


def min_int_dtype(max_abs_value: int | float, min_value: int | float, allow_uints: bool = True) -> np.dtype:
    if allow_uints and min_value > 0:
        return min_uint_dtype(max_abs_value)

    # Find the smallest dtype that can represent all values in the array
    for dtype in DTYPE_INT:
        if max_abs_value <= np.iinfo(dtype).max:
            return dtype

    raise ValueError(f"Number is larger than '{DTYPE_INT[-1]}' and can't be represented.")


def min_int_dtype_array(array: np.ndarray, allow_uints: bool = True) -> np.dtype:
    return min_int_dtype(max(np.max(array), np.abs(np.min(array))), np.min(array), allow_uints)


def min_float_dtype(max_abs_value: float) -> np.dtype:
    # Find the smallest dtype that can represent all values in the array
    for dtype in DTYPE_FLOAT:
        if max_abs_value <= np.iinfo(dtype).max:
            return dtype

    raise ValueError(f"Number is larger than '{DTYPE_FLOAT[-1]}' and can't be represented.")


def min_float_dtype_array(array: np.ndarray) -> np.dtype:
    return min_float_dtype(max(np.max(array), np.abs(np.min(array))))


def set_minimum_dtype(array: np.ndarray) -> np.ndarray:
    if array.dtype in DTYPE_UINT or array.dtype in DTYPE_INT:
        dtype_ = min_int_dtype_array(array)
        return array.astype(dtype_)

    if array.dtype in DTYPE_FLOAT:
        dtype_ = min_float_dtype_array(array)
        return array.astype(dtype_)

    raise ValueError("Not supported dtype.")


def map_discrete_x_axis(x_new: np.ndarray, x_old: np.ndarray, y_old: np.ndarray, ignore_issues: bool = False) \
        -> np.ndarray:
    y_new = np.zeros_like(x_new, dtype=y_old.dtype)
    index = get_index_of_values_in_common(x_new, x_old)
    if ignore_issues and any(i is None for i in index):
        zip_list = ((y, i) for y, i in zip(y_old, index) if i is not None)
        y_old, index = zip(*zip_list)
    else:
        if any(i is None for i in index):
            raise ValueError("Mapping not.")
    y_new[list(index)] = y_old

    return y_new


def map_discrete_x_axis_2D(x_new: np.ndarray, x_old: np.ndarray, y_old: np.ndarray, ignore_issues: bool = False) \
        -> np.ndarray:
    y_new = np.zeros((y_old.shape[0], len(x_new)), dtype=y_old.dtype)
    index = get_index_of_values_in_common(x_new, x_old)

    for i, index_ in enumerate(index):
        if ignore_issues and index_ is None:
            continue
        else:
            if index_ is None:
                raise ValueError("Mapping not.")
        y_new[:, index_] = y_old[:, i]

    return y_new


def get_all_unique_values_within_tolerance(arrays: Sequence[np.ndarray], **kwargs) -> np.ndarray:
    if len(kwargs) == 0:
        kwargs = {"rtol": 0.01}
    combined_array = np.concatenate(arrays)
    sorted_array = np.sort(combined_array)

    # Compare each value to the last added unique value
    unique_values = [sorted_array[0]]
    for value in sorted_array[1:]:
        if not np.isclose(value, unique_values[-1], **kwargs):
            unique_values.append(value)

    return np.array(unique_values)


def rescale_array(
        arr: np.ndarray,
        new_lb: int | float,
        new_up: int | float,
        old_lb: int | float = 0,
        old_up: int | float = 1
) -> np.ndarray:
    return new_lb + (new_up - new_lb) * (arr - old_lb) / (old_up - old_lb)


def get_slope(y: np.ndarray) -> float:
    """Compute the slope of the last few values in `norms` using linear regression."""
    if len(y) < 2:
        return 0  # Not enough points to compute a meaningful slope

    x = np.arange(len(y))
    A = np.vstack([x, np.ones_like(x)]).T
    slope, _ = np.linalg.lstsq(A, y, rcond=None)[0]

    return slope
