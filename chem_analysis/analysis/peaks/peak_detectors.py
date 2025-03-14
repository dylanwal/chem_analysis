from typing import Any

import numpy as np
from scipy.signal import find_peaks as scipy_find_peaks
from scipy.ndimage import grey_dilation, grey_opening

from chem_analysis.processing.baseline.morphological import estimate_window
from chem_analysis.analysis.peaks.base_classes import Detector
import chem_analysis.utils.math as math_utils


class ScipyPeakFinder(Detector):
    def __init__(self,
                 height: float | np.ndarray | None = None,
                 threshold: float | np.ndarray | None = None,
                 distance: float | None = None,
                 prominence: float | np.ndarray | None = None,
                 width: float | np.ndarray | None = None,
                 wlen: int | None = None,
                 rel_height: float | None = None,
                 plateau_size: float | np.ndarray | None = None,
                 ):
        self.height = height
        self.threshold = threshold
        self.distance = distance
        self.prominence = prominence
        self.width = width
        self.wlen = wlen
        self.rel_height = rel_height
        self.plateau_size = plateau_size

    def get_kwargs(self) -> dict[str, Any]:
        keys = ['height', 'threshold', 'distance', 'prominence', 'width', 'wlen', 'rel_height', 'plateau_size']
        return {arg: getattr(self, arg) for arg in keys if getattr(self, arg) is not None}

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        index, _ = scipy_find_peaks(y, **self.get_kwargs())
        return index


class MaxValues(Detector):
    def __init__(self, n: int = 1):
        """

        Parameters
        ----------
        n: number of max values

        """
        self.n = n

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        index = np.argsort(y)
        return index[:self.n]


class LocalMax(Detector):
    def __init__(self):
        pass

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        from scipy.signal._peak_finding_utils import _local_maxima_1d
        index, _, _ = _local_maxima_1d(y)
        return index


class Derivative(Detector):
    def __init__(self):
        pass

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Find peaks using the first and second derivatives."""
        dy = np.gradient(y, x)  # First derivative
        ddy = np.gradient(dy, x)  # Second derivative

        zero_crossings = np.where(np.diff(np.sign(dy)) < 0)[0]  # First derivative crosses zero downward
        peaks = [i for i in zero_crossings if ddy[i] < 0]

        return np.array(peaks, dtype=int)  # Peak x and y values


def estimate_scales(y: np.ndarray, window: int | None = None, min_scale: int = 1, max_scale: int =50):
    """Estimate optimal scales based on peak width."""
    if window is None:
        window = int(estimate_window(y) / 5)
    y_ = grey_opening(y, window)
    bounds = math_utils.find_consecutive_regions(y_)
    if bounds is None:
        return np.arange(min_scale, max_scale)

    median_width = np.median(np.diff(bounds, axis=1).flatten())
    return np.arange(min_scale, int(min(max_scale, median_width * 2)))


# import pywt
#
#
# class Wavelet(Detector):
#     def __init__(self,
#                  wavelet: str = 'mexh',
#                  threshold: float = 0.1,
#                  scales: np.ndarray | None = None,
#                  ):
#         self.wavelet = wavelet
#         self.threshold = abs(threshold)
#         self.scales = scales
#         self._bounds = None
#
#     def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
#         if self.scales is None:
#             # self.scales = estimate_scales(y)
#             min_freq = 20/len(y)
#             max_freq = len(y)/20/len(y)
#             self.scales = np.linspace(1, 64, 32)
#         # Compute the continuous wavelet transform (CWT) using PyWavelets
#         coeffs, _ = pywt.cwt(y, self.scales, self.wavelet)
#
#         # Use the scale that gives the strongest response
#         cwt_sum = np.sum(np.abs(coeffs), axis=0)
#         peaks = Derivative().run_xy(x, cwt_sum)
#
#         index = np.argsort(y[peaks])
#         peak_apex = peaks[index[:]]
#
#         # Find the beginning (left boundary) and end (right boundary)
#         for peak in peak_apex:
#             left_idx = np.nonzero(y[:peak_apex] < y[peak_apex] * self.threshold)[0]
#             right_idx = np.nonzero(y[peak_apex:] < y[peak_apex] * self.threshold)[0]
#
#             peak_start = left_idx[-1] if len(left_idx) > 0 else 0
#             peak_end = peak_apex + right_idx[0] if len(right_idx) > 0 else len(y) - 1
#             self._bounds = np.column_stack((peak_start, peak_end))
#
#         return peak_apex


class Morphological(Detector):
    def __init__(self,
                 min_: int | float = None,
                 max_: int | float = None,
                 window: int | None = None,
                 mode: str = "index",
                 auto_div: int | float | None = None,
                 ):
        """

        Parameters
        ----------
        min_
        max_
        window
        mode:
            "span"  x distance
            "index" index
        auto_div:
            divisor for auto window algorithm
        """
        if not (mode == 'span' or mode == 'index'):
            raise ValueError("Type must be 'span' or 'index'")
        self.min_ = min_
        self.max_ = max_
        self.window = window
        self.mode = mode
        self.auto_div = auto_div or 5
        self._bounds = None

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window is None:
            self.window = int(estimate_window(y) / self.auto_div)

        y_ = grey_dilation(y, 80)
        bounds = math_utils.find_consecutive_regions(y_)

        # compute plateau that are big or small enough
        if self.mode == "index":
            span = np.diff(bounds, axis=1).flatten()
        elif self.mode == "span":
            span = np.diff(x[bounds], axis=1).flatten()
        else:
            raise ValueError("Mode must be 'span' or 'index'")  # checked earlier

        mask = np.ones_like(span, dtype=bool)
        if self.min_ is not None:
            mask &= span >= self.min_
        if self.max_ is not None:
            mask &= span <= self.max_
        bounds = bounds[mask]
        self._bounds = bounds

        # get middle index
        index = math_utils.get_middle_index(self._bounds)

        # remove 0 thin peaks
        index = index[self._bounds[:, 0] != index]
        return index


class MovingAverage(Detector):
    def __init__(self,
                 lag: int = 5,
                 threshold: float | float = 3,
                 influence: float | int = 0.05,
                 mode: str = "positive",
                 ):
        """

        Parameters
        ----------
        lag
        threshold
        influence
        mode:
            "positive": only positive peaks
            "negative": only negative peaks
            "both": both positive and negative peaks
        """
        self.lag = lag
        self.threshold = threshold
        self.influence = influence
        self._bounds = None
        self.mode = mode

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        lag, threshold, influence = self.lag, self.threshold, self.influence

        filteredY = y.copy()
        avgFilter = np.empty_like(y)
        stdFilter = np.empty_like(y)
        signals = np.zeros_like(y)

        avgFilter[lag - 1] = np.mean(y[:lag])
        stdFilter[lag - 1] = np.std(y[:lag])

        for i in range(lag, len(y)):
            if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter[i-1]:
                signals[i] = np.sign(y[i] - avgFilter[i-1])  # 1 if above, -1 if below
                filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            else:
                signals[i] = 0
                filteredY[i] = y[i]

            # Update rolling statistics efficiently
            window = filteredY[(i - lag + 1):(i + 1)]
            avgFilter[i] = np.mean(window)
            stdFilter[i] = np.std(window)

        # plotting
        # import plotly.graph_objects as go
        # fig = go.Figure()
        # fig.add_scatter(x=x, y=y, name='raw')
        # fig.add_scatter(x=x, y=avgFilter - stdFilter*threshold, name='lowwer')
        # fig.add_scatter(x=x, y=avgFilter + stdFilter*threshold, fill='tonexty', name='upper')
        # fig.add_scatter(x=x, y=filteredY, name='filteredY')
        # fig.add_scatter(x=x, y=signals*50_000, name='sig')
        # fig.show()

        signals[-1] = 0  # ensure it ends down

        if self.mode == "positive":
            signals[signals == -1] = 0
        elif self.mode == "negative":
            signals[signals == 1] = 0
            signals[signals == -1] = 1
        elif self.mode == "both":
            signals[signals == -1] = 1
        else:
            raise ValueError("Mode must be 'positive' or 'negative' or 'both'")

        up = np.nonzero(np.diff(signals) > 0)[0]
        downs = np.nonzero(np.diff(signals) < 0)[0]
        self._bounds = np.column_stack((up, downs))
        index = math_utils.get_middle_index(self._bounds)

        # remove 0 thin peaks
        index = index[self._bounds[:, 0] != index]
        return index
