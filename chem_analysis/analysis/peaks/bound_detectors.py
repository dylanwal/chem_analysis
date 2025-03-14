import numpy as np
from scipy.ndimage import grey_dilation, grey_opening

import chem_analysis.utils.math as math_utils
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.processing.processing_method import Smoothing
from chem_analysis.processing.baseline.morphological import estimate_window


def _apply_smoother(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray,
        smoother: Smoothing,
        adjust_index: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if smoother is not None:
        x, y = np.copy(x), np.copy(y)
        x, y = smoother.run(x, y)

        # smoothing can cause previous max not to be the peak anymore; so move it
        if adjust_index:
            index = np.copy(index)
            for i, index_ in enumerate(index):
                index[i] = math_utils.find_local_min(-1 * y, index_)
            index = np.unique(index)

    return x, y, index


def bounds_max_slope(
        signal: Signal,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,
) -> np.ndarray:
    return bounds_max_slope_xy(signal.x, signal.y, index, smoother, adjust_index)


def bounds_max_slope_xy(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,
) -> np.ndarray:
    if smoother is not None:
        x, y, index = _apply_smoother(x, y, index, smoother, adjust_index)

    bounds = np.empty((len(index), 2), dtype=int)
    for i, idx in enumerate(index):
        bounds[i] = math_utils.find_max_slope(y, idx, len(y))
    return bounds


def bounds_first_increase(
        signal: Signal,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,
) -> np.ndarray:
    return bounds_first_increase_xy(signal.x, signal.y, index, smoother, adjust_index)


def bounds_first_increase_xy(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,

) -> np.ndarray:
    if smoother is not None:
        x, y, index = _apply_smoother(x, y, index, smoother, adjust_index)

    bounds = np.empty((len(index), 2), dtype=int)
    for i, idx in enumerate(index):
        bounds[i] = math_utils.find_first_increase(y, idx, len(y))
    return bounds


def bounds_first_increase_zero(
        signal: Signal,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,
) -> np.ndarray:
    return bounds_first_increase_zero_xy(signal.x, signal.y, index, smoother, adjust_index)


def bounds_first_increase_zero_xy(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray,
        smoother: Smoothing = None,
        adjust_index: bool = True,
) -> np.ndarray:
    if smoother is not None:
        x, y, index = _apply_smoother(x, y, index, smoother, adjust_index)

    bounds = np.empty((len(index), 2), dtype=int)
    for i, idx in enumerate(index):
        bounds[i] = math_utils.find_first_increase_zero(y, idx, len(y))
    return bounds


def bounds_morph_dilation(
        signal: Signal,
        index: np.ndarray | None = None,
        window: int | None = None,
        mode: str = "plateau",
        auto_div: int | float | None = None,
) -> np.ndarray:
    return bounds_morph_dilation_xy(signal.x, signal.y, index, window, mode, auto_div)


def bounds_morph_dilation_xy(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray | None = None,
        window: int | None = None,
        mode: str = "plateau",
        auto_div: int | float | None = None,
) -> np.ndarray:
    """

    Parameters
    ----------
    x:
    y:
    window
    index:
        if you want the bounds for specific peaks; provide the index of peak max
    mode:
        "plateau"
        "max_slope"
        "first_increase"
        "first_increase_zero"
    auto_div

    Returns
    -------
    bounds
    [[left, right], [left, right], ...]

    """
    if window is None:
        if auto_div is None:
            auto_div = 2
        window = int(estimate_window(y) / auto_div)

    y_ = grey_dilation(y, window)
    if mode == "plateau":
        if index is None:
            bounds = math_utils.find_consecutive_regions(y_)
        else:
            bounds = np.empty((len(index), 2), dtype=int)
            for i, idx in enumerate(index):
                bounds[i] = math_utils.find_first_change_in_value(y_, idx, len(y))
    elif mode == "max_slope":
        if index is None:
            index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_max_slope(y_, idx, len(y))
    elif mode == "first_increase":
        if index is None:
            index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_first_increase(y_, idx, len(y))
    elif mode == "first_increase_zero":
        if index is None:
            index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_first_increase_zero(y_, idx, len(y), exclude=True)

    else:
        raise ValueError("Not supported 'mode'.")

    return bounds


def bounds_morph_opening(
        signal: Signal,
        index: np.ndarray | None = None,
        window: int | None = None,
        auto_div: int | float | None = None,
) -> np.ndarray:
    return bounds_morph_opening(signal.x, signal.y, index, window, auto_div)


def bounds_morph_opening_xy(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray | None = None,
        window: int | None = None,
        auto_div: int | float | None = None,
) -> np.ndarray:
    """

    Parameters
    ----------
    x:
    y:
    window
    index:
        if you want the bounds for specific peaks; provide the index of peak max
    auto_div

    Returns
    -------
    bounds
    [[left, right], [left, right], ...]

    """
    if window is None:
        if auto_div is None:
            auto_div = 3
        window = int(estimate_window(y) / auto_div)

    y_ = grey_opening(y, window)
    if index is None:
        bounds = math_utils.find_consecutive_regions(y_)
    else:
        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_first_change_in_value(y_, idx, len(y))
    return bounds


# def bounds_wavelet(
#         signal: Signal,
#         index: np.ndarray | None = None,
#         wavelet: str = 'mexh',
#         threshold: float = 0.1,
#         scales: np.ndarray | None = None,
# ) -> np.ndarray:
#     return bounds_wavelet_xy(signal.x, signal.y, index, wavelet, threshold, scales)
#
#
# def bounds_wavelet_xy(
#         x: np.ndarray,
#         y: np.ndarray,
#         index: np.ndarray | None = None,
#         wavelet: str = 'mexh',
#         threshold: float = 0.1,
#         scales: np.ndarray | None = None,
# ) -> np.ndarray:
#     """
#
#     Parameters
#     ----------
#     x:
#     y:
#
#
#     Returns
#     -------
#     bounds
#     [[left, right], [left, right], ...]
#
#     """
#     method = Wavelet(wavelet, threshold, scales)
#     method.run_xy(x, y)
#     bounds_ = method._bounds()
#
#     if index is None:
#         return bounds_
#
#     index_ = np.empty(len(index), dtype=int)
#     for i, idx in enumerate(index):
#         mask = (bounds_[0] < idx) & (idx < bounds_[1])
#         indexes = np.nonzero(mask)[0]
#         if len(indexes) < 0:
#             index_[i] = -1
#         else:
#             index_[i] = indexes[0]
#
#     bounds = bounds_[index_]
#     bounds[index_ == -1] = [0, 0]
#
#     return bounds
