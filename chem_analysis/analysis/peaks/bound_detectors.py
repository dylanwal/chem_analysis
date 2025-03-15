import numpy as np
from scipy.ndimage import grey_dilation, grey_opening

import chem_analysis.utils.math as math_utils
from chem_analysis.processing.processing_method import Smoothing
from chem_analysis.processing.baseline.morphological import estimate_window
from chem_analysis.analysis.peaks.base_classes import BoundDetector


def _apply_smoother(
        x: np.ndarray,
        y: np.ndarray,
        index: np.ndarray,
        smoother: Smoothing,
        adjust_index: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if smoother is not None:
        x, y = np.copy(x), np.copy(y)
        x, y = smoother.run_xy(x, y)

        # smoothing can cause previous max not to be the peak anymore; so move it
        if adjust_index:
            index = np.copy(index)
            for i, index_ in enumerate(index):
                index[i] = math_utils.find_local_min(-1 * y, index_)
            index = np.unique(index)

    return x, y, index


class BoundMaxSlope(BoundDetector):
    def __init__(self,
                 adjust_index: bool = True,
                 smoother: Smoothing = None
                 ):
        self.smoother = smoother
        self.adjust_index = adjust_index

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.smoother is not None:
            x, y, index = _apply_smoother(x, y, index, self.smoother, self.adjust_index)

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_max_slope(y, idx, len(y))
        return bounds


class BoundFirstIncrease(BoundDetector):
    def __init__(self,
                 adjust_index: bool = True,
                 smoother: Smoothing = None
                 ):
        self.smoother = smoother
        self.adjust_index = adjust_index

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.smoother is not None:
            x, y, index = _apply_smoother(x, y, index, self.smoother, self.adjust_index)

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_first_increase(y, idx, len(y))
        return bounds


class BoundFirstIncreaseZero(BoundDetector):
    def __init__(self,
                 adjust_index: bool = True,
                 smoother: Smoothing = None
                 ):
        self.smoother = smoother
        self.adjust_index = adjust_index

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.smoother is not None:
            x, y, index = _apply_smoother(x, y, index, self.smoother, self.adjust_index)

        bounds = np.empty((len(index), 2), dtype=int)
        for i, idx in enumerate(index):
            bounds[i] = math_utils.find_first_increase_zero(y, idx, len(y))
        return bounds


class BoundMorphDilation(BoundDetector):
    def __init__(self,
                 window: int | None = None,
                 mode: str = "plateau",
                 auto_div: int | float | None = None,
                 adjust_index: bool = True,
                 smoother: Smoothing = None
                 ):
        """

        Parameters
        ----------
        window

        mode:
            "plateau"
            "max_slope"
            "first_increase"
            "first_increase_zero"
        auto_div
        """
        self.window = window
        self.mode = mode
        self.auto_div = auto_div
        self.smoother = smoother
        self.adjust_index = adjust_index

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.smoother is not None:
            x, y, index = _apply_smoother(x, y, index, self.smoother, self.adjust_index)

        if self.window is None:
            if self.auto_div is None:
                self.auto_div = 2
            self.window = int(estimate_window(y) / self.auto_div)

        y_ = grey_dilation(y, self.window)
        if self.mode == "plateau":
            if index is None:
                bounds = math_utils.find_consecutive_regions(y_)
            else:
                bounds = np.empty((len(index), 2), dtype=int)
                for i, idx in enumerate(index):
                    bounds[i] = math_utils.find_first_change_in_value(y_, idx, len(y))
        elif self.mode == "max_slope":
            if index is None:
                index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

            bounds = np.empty((len(index), 2), dtype=int)
            for i, idx in enumerate(index):
                bounds[i] = math_utils.find_max_slope(y_, idx, len(y))
        elif self.mode == "first_increase":
            if index is None:
                index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

            bounds = np.empty((len(index), 2), dtype=int)
            for i, idx in enumerate(index):
                bounds[i] = math_utils.find_first_increase(y_, idx, len(y))
        elif self.mode == "first_increase_zero":
            if index is None:
                index = math_utils.get_middle_index(math_utils.find_consecutive_regions(y_))

            bounds = np.empty((len(index), 2), dtype=int)
            for i, idx in enumerate(index):
                bounds[i] = math_utils.find_first_increase_zero(y_, idx, len(y), exclude=True)

        else:
            raise ValueError("Not supported 'mode'.")

        return bounds


class BoundMorphOpening(BoundDetector):
    def __init__(self,
                 window: int | None = None,
                 auto_div: int | float | None = None,
                 adjust_index: bool = True,
                 smoother: Smoothing = None
                 ):
        self.window = window
        self.auto_div = auto_div
        self.smoother = smoother
        self.adjust_index = adjust_index

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.smoother is not None:
            x, y, index = _apply_smoother(x, y, index, self.smoother, self.adjust_index)

        if self.window is None:
            if self.auto_div is None:
                self.auto_div = 3
            self.window = int(estimate_window(y) / self.auto_div)

        y_ = grey_opening(y, self.window)
        if index is None:
            bounds = math_utils.find_consecutive_regions(y_)
        else:
            bounds = np.empty((len(index), 2), dtype=int)
            for i, idx in enumerate(index):
                bounds[i] = math_utils.find_first_change_in_value(y_, idx, len(y))
        return bounds



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
