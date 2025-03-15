from typing import Iterable, Sequence
import numpy as np

from chem_analysis.analysis.peaks.base_classes import PeakFilter
import chem_analysis.utils.math as math_utils


def _slice_to_mask(slice_: slice, index: np.ndarray, max_index: int) -> np.ndarray:
    if slice_.start is None:
        start = 0
    else:
        start = slice_.start
    if slice_.stop is None:
        stop = max_index
    else:
        stop = slice_.stop

    sub_mask = index > start
    sub_mask &= index < stop
    return sub_mask


class FilterSlices(PeakFilter):
    def __init__(self, slices: slice | Iterable[slice], invert: bool = False):
        """

        Parameters
        ----------
        slices:
            valid format slice(start, stop) or (slice(start, stop), slice(start, stop), ...)
        invert:
            False: slices are valid locations for peaks
            True: slices are invalid locations for peaks
        """
        if not isinstance(slices[0], Iterable):
            slices = [slices]

        self.slices = slices
        self.invert = invert

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        mask = np.zeros_like(index, dtype=bool)
        for slice_ in self.slices:
            mask = np.logical_or(mask, _slice_to_mask(slice_, index, len(y)))

        if self.invert:
            mask = np.logical_not(mask)
        return index[mask]


class FilterSpans(PeakFilter):
    def __init__(self,
                 spans: Sequence[float | None] | Iterable[Sequence[float | None]],  # Sequence of length 2
                 invert: bool = False
                 ):
        """

        Parameters
        ----------
        spans:
            valid format (start, stop) or ((start, stop), (start, stop), ...)
        invert:
            False: slices are valid locations for peaks
            True: slices are invalid locations for peaks
        """
        if not isinstance(spans[0], Iterable):
            spans = [spans]

        self.spans = spans
        self.invert = invert

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        x_ = x[index]
        mask = np.ones_like(index, dtype=bool)
        for x_span in self.spans:
            left = x_span[0] or x[0]
            right = x_span[1] or x[-1]
            sub_mask = x_ > left
            sub_mask &= x_ < right
            mask &= sub_mask

        if self.invert:
            mask = np.logical_not(mask)
        return index[mask]


class FilterHeight(PeakFilter):
    def __init__(self,
                 min_abs: float = None,
                 max_abs: float = None,
                 min_rel: float = None,
                 max_rel: float = None,
                 ):
        """

        Parameters
        ----------
        min_abs:
            with respect to y
        max_abs
            with respect to y
        min_rel
            with respect to max(y)
            typically between [0, 1)
        max_rel:
            with respect to max(y)
            typically between [0, 1) # with 1=tallest peak
        """
        if min_abs is None and max_abs is None and min_rel is None and max_rel is None:
            raise ValueError("At least on height parameter is required")
        self.min_abs = min_abs
        self.max_abs = max_abs
        self.min_rel = min_rel
        self.max_rel = max_rel

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        heights = y[index]

        mask = np.ones_like(index, dtype=bool)
        if self.min_abs is not None:
            np.logical_and(heights >= self.min_abs, mask, out=mask)
        if self.max_abs is not None:
            np.logical_and(heights <= self.max_abs, mask, out=mask)
        if self.min_rel is not None:
            np.logical_and(heights >= self.min_rel * np.max(y), mask, out=mask)
        if self.max_rel is not None:
            np.logical_and(heights <= self.max_rel * np.max(y), mask, out=mask)

        return index[mask]


class FilterHeightLocal(PeakFilter):
    def __init__(self,
                 min_abs: float = None,
                 max_abs: float = None,
                 min_rel: float = None,
                 max_rel: float = None,
                 mode: str = 'peak',
                 window: int | float = 1,
                 ):
        """

        Parameters
        ----------
        min_abs
        max_abs
        min_rel
        max_rel
        mode:
            "peak": window is peak index
            "x": window is x range
        window:
            will look left and right window size
        """
        if min_abs is None and max_abs is None and min_rel is None and max_rel is None:
            raise ValueError("At least on height parameter is required")
        if not (mode == 'peak' or mode == 'x'):
            raise ValueError("Type must be 'peak' or 'x'")
        self.min_abs = min_abs
        self.max_abs = max_abs
        self.min_rel = min_rel
        self.max_rel = max_rel
        self.mode = mode
        self.window = window

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        y_ = y[index]
        x_ = x[index]

        mask = np.ones_like(index, dtype=bool)
        for i in range(len(index)):
            if self.mode == 'peak':
                sub_mask = np.zeros_like(index, dtype=bool)
                sub_mask[max(0, i-self.window):i] = True
                sub_mask[i:min(len(index), i + self.window)+1] = True
            elif self.mode == 'x':
                sub_mask = x_ > x_[i] - self.window
                sub_mask &= x_ < x_[i] + self.window
            else:
                continue

            sub_mask[i] = False  # remove the peak itself
            if not np.any(sub_mask):
                continue

            distance = np.abs(y_[sub_mask] - y_[i])
            if self.min_abs is not None:
                mask[i] &= np.all(distance >= self.min_abs)
            if self.max_abs is not None:
                mask[i] &= np.all(distance <= self.max_abs)
            if self.min_rel is not None:
                mask[i] &= np.all(distance >= self.min_rel * np.max(y))
            if self.max_rel is not None:
                mask[i] &= np.all(distance <= self.max_rel * np.max(y))

        return index[mask]


class FilterSpacing(PeakFilter):
    def __init__(self,
                 min_: float = None,
                 max_: float = None  # not sure if we need max??
                 ):
        """
        horizontal spacing
        smallest y height values removed first
        Parameters
        ----------
        min_
        max_
        """
        if min_ is None and max_ is None:
            raise ValueError("At least one spacing parameter is required")
        self.min_ = min_
        self.max_ = max_

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        sorted_index = np.argsort(y[index])  # sort by height
        mask = np.ones_like(sorted_index, dtype=bool)
        for sort_index in sorted_index:
            x_peak = x[index[sort_index]]

            if sort_index != 0:
                for i in range(sort_index - 1, 0, -1):
                    if mask[i]:
                        x_left = x[index[i]]
                        if self.min_ is not None and abs(x_peak - x_left) < self.min_:
                            mask[sort_index] = False
                        if self.max_ is not None and abs(x_peak - x_left) > self.max_:
                            mask[sort_index] = False
                        break

            if sort_index < len(sorted_index) - 1:
                for i in range(sort_index + 1, len(sorted_index)):
                    if mask[i]:
                        x_right = x[index[i]]
                        if self.min_ is not None and abs(x_peak - x_right) < self.min_:
                            mask[sort_index] = False
                        if self.max_ is not None and abs(x_peak - x_right) > self.max_:
                            mask[sort_index] = False

        return index[mask]


class FilterProminence(PeakFilter):
    def __init__(self,
                 min_: float = None,
                 max_: float = None,
                 window: int = -1,
                 ):
        """
        horizontal spacing
        smallest y height values removed first
        Parameters
        ----------
        min_
        max_
        window
        """
        if min_ is None and max_ is None:
            raise ValueError("At least one spacing parameter is required")
        self.min_ = min_
        self.max_ = max_
        self.window = window

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        from scipy.signal._peak_finding_utils import _peak_prominences
        result = _peak_prominences(y, np.intp(index), np.intp(self.window))

        prominences = result[0]
        mask = np.ones_like(index, dtype=bool)
        if self.min_ is not None:
            mask &= (self.min_ <= prominences)
        if self.max_ is not None:
            mask &= (prominences <= self.max_)

        return index[mask]


class FilterWidth(PeakFilter):
    def __init__(self,
                 min_: float = None,
                 max_: float = None,
                 height: float = 0.5,
                 mode: str = 'span',
                 ):
        """
        peak width filter

        Parameters
        ----------
        min_
        max_
        height:
            relative to max height that it will try to compute width.
            may be less if height increase before (0.5*peak height) reached
        mode:
            "span"  x distance
            "index" index
        """
        if min_ is None and max_ is None:
            raise ValueError("At least one spacing parameter is required")
        if not (mode == 'span' or mode == 'index'):
            raise ValueError("Type must be 'span' or 'index'")
        self.min_ = min_
        self.max_ = max_
        self.height = height
        self.mode = mode

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        mask = np.ones_like(index, dtype=bool)
        for i, idx in enumerate(index):
            idx = int(idx)
            cut_off_height = y[idx] * self.height

            # left
            for ii in range(idx - 1, 0, -1):
                if y[ii] < y[ii - 1] or y[ii] < cut_off_height:  # slope increase or height reached
                    left = ii
                    break
            else:
                left = 0

            # right
            for ii in range(idx + 1, len(y) + 1):
                if y[ii] < y[ii + 1] or y[ii] < cut_off_height:  # slope increase or height reached
                    right = ii
                    break
            else:
                right = len(y)

            if self.mode == 'span':
                span = x[right] - x[left]
            elif self.mode == 'index':
                span = right - left
            else:
                raise ValueError("Type must be 'span' or 'index'")

            if self.min_ is not None and span < self.min_:
                mask[i] &= False
            if self.max_ is not None and span > self.max_:
                mask[i] &= False

        return index[mask]


from chem_analysis.processing.baseline.morphological import estimate_window
from scipy.ndimage import grey_opening


class FilterWidthMorphological(PeakFilter):
    def __init__(self,
                 min_: float = None,
                 max_: float = None,
                 window: int | None = None,
                 mode: str = 'span',
                 auto_div: int | float | None = None,
                 ):
        """
        horizontal spacing
        smallest y height values removed first
        Parameters
        ----------
        min_
        max_
        window:
            window used in morphological filter
        mode:
            "span"  x distance
            "index" index
        auto_div:
            divisor to tune the auto window algorithm
        """
        if min_ is None and max_ is None:
            raise ValueError("At least one spacing parameter is required")
        if not (mode == 'span' or mode == 'index'):
            raise ValueError("Type must be 'span' or 'index'")
        self.min_ = min_
        self.max_ = max_
        self.window = window
        self.mode = mode
        self.auto_div = auto_div or 5

    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        if self.window is None:
            self.window = int(estimate_window(y)/self.auto_div)

        y_ = grey_opening(y, self.window)
        # all peaks will now have flat tops; now find width of flat tops

        mask = np.ones_like(index, dtype=bool)
        for i, idx in enumerate(index):
            left, right = math_utils.find_first_change_in_value(y_, int(idx), len(y))

            if self.mode == 'span':
                span = x[right] - x[left]
            elif self.mode == 'index':
                span = right - left
            else:
                raise ValueError("Type must be 'span' or 'index'")

            if self.min_ is not None and span < self.min_:
                mask[i] &= False
            if self.max_ is not None and span > self.max_:
                mask[i] &= False

        return index[mask]
