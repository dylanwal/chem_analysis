from typing import Union

import numpy as np


def rolling_value(x: np.ndarray, peak_index: int = None, peak_value: Union[float, int] = None,
                  sensitivity: float = 0.03, cut_off: Union[None, float, int] = 0) -> tuple[int, int]:
    """

    Parameters
    ----------
    x: np.ndarray
        1D data to find bounds
    peak_index: int
        index of peak max
    peak_value: int, float
        value of peak max
    sensitivity: float
        How much it can go up before triggering a bound detection, fraction of max height
    cut_off: Union[None, float, int]
        When to stop if never goes to zero, fraction of max height

    Returns
    -------
    lb: int
        index of lower bound
    up: int
        index of upper bound
    """
    if peak_index is None:
        raise ValueError("Peak index or peak value are required.")
    if peak_value:
        peak_index_ = np.where(x == peak_value)
        if len(peak_index_[0]) == 0:
            raise ValueError(f"Peak not found in data. (peak_value: {peak_value})")
        elif len(peak_index_[0]) >= 2:
            raise ValueError(f"Multiple peaks have this peak_value. (peak_value: {peak_value})")

        peak_index = peak_index_[0]
    else:
        peak_value = x[peak_index]

    # set up cutoffs
    if cut_off is not None:
        cut_off_value = peak_value * cut_off
    else:
        cut_off_value = -10e-50
    if cut_off is not None:
        sensitivity_value = peak_value * sensitivity
    else:
        sensitivity_value = 10e-50

    # lower bound
    if peak_index == 0:
        lb = peak_index
    elif peak_index < 5:
        lb = np.argmin(x[:peak_index])
    else:
        min_ = x[peak_index]
        min_index = peak_index
        for i, v in enumerate(np.flip(x[:peak_index])):
            if v < min_:
                min_ = v
                min_index = peak_index - i
                if min_ < cut_off_value:
                    break
            if v - x[min_index] > sensitivity_value:
                break

        lb = min_index-1

    # upper bound
    if peak_index == len(x):
        ub = peak_index
    elif len(x) - peak_index < 5:
        ub = np.argmin(x[peak_index:])
    else:
        min_ = x[peak_index]
        min_index = peak_index
        for i, v in enumerate(x[peak_index:]):
            if v < min_:
                min_ = v
                min_index = peak_index + i
                if min_ < cut_off_value:
                    break
            if v - x[min_index] > sensitivity_value:
                break

        ub = min_index

    return lb, ub


bound_detection_methods = {
    "rolling_value": rolling_value
}
