from typing import Any

import numpy as np
from scipy.signal import find_peaks as scipy_find_peaks

from chem_analysis.analysis.peak_picking.base_classes import Detector


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
        return {arg:getattr(self, arg) for arg in keys if getattr(self, arg) is not None}

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


def estimate_scales(signal, min_scale=1, max_scale=50):
    """Estimate optimal scales based on peak width."""
    peaks, properties = scipy_find_peaks(signal, prominence=0.1)
    if len(peaks) > 0:
        peak_widths = properties["prominences"]  # Use prominence as a proxy for width
        median_width = np.median(peak_widths) if len(peak_widths) > 0 else 10
        scales = np.arange(min_scale, int(min(max_scale, median_width * 2)))
    else:
        scales = np.arange(min_scale, max_scale)  # Default if no peaks found
    return scales


class Wavelet(Detector):
    def __init__(self):
        wavelet='mexh'
        pass

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        scales = estimate_scales(y)
        # Compute the continuous wavelet transform (CWT) using PyWavelets
        coeffs, _ = pywt.cwt(y, scales, wavelet)

        # Use the scale that gives the strongest response
        cwt_sum = np.sum(np.abs(coeffs), axis=0)
        peaks, _ = find_peaks(cwt_sum)

        if len(peaks) == 0:
            return None, None, None

        index = np.argsort(y[peaks])
        peak_apex = peaks[index[:]]
        return peak_apex

        # # Find the beginning (left boundary) and end (right boundary)
        # left_idx = np.where(y[:peak_apex] < y[peak_apex] * 0.1)[0]
        # right_idx = np.where(y[peak_apex:] < y[peak_apex] * 0.1)[0]
        #
        # peak_start = left_idx[-1] if len(left_idx) > 0 else 0
        # peak_end = peak_apex + right_idx[0] if len(right_idx) > 0 else len(y) - 1
        #
        # return peak_start, peak_apex, peak_end



