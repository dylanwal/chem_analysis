import numpy as np

from chem_analysis.processing.smoothing.savitzky_golay import Smoothing


def detect_outliers(data: np.ndarray, m: float):
    dist_from_median = np.abs(data - np.median(data))
    median_deviation = np.median(dist_from_median)
    if median_deviation != 0:
        scale_distances_from_median = dist_from_median / median_deviation
        return scale_distances_from_median > m  # True is an outlier

    return np.zeros_like(data)  # no outliers


def window_calc(data: np.ndarray, pos: int, m: float) -> float:
    outlier = detect_outliers(data, m)
    if outlier[pos]:
        return float(np.median(data[np.invert(outlier)]))
    else:
        return data[pos]


def rolling_window_main(window: float, m: float, y):
    out = np.empty_like(y)

    if window % 2 != 0:
        window += 1

    span = int(window / 2)
    for i in range(len(y)):
        if i < span:  # left edge
            out[i] = window_calc(y[:window], i, m)
        elif i > len(y) - span:  # right edge
            out[i] = window_calc(y[window:], int(i - (len(y) - window)), m)
        else:  # middle
            out[i] = window_calc(y[i - span:i + span], span, m)

    return out


class RollingWindow(Smoothing):
    def __init__(self, window: int = 20, m: float = 2, temporal_processing: int = 1):
        super().__init__(temporal_processing)
        self.window = window
        self.m = m

    def run_xy(self, x: np.ndarray, y: np.ndarray, ) -> tuple[np.ndarray, np.ndarray]:
        return x, rolling_window_main(self.window, self.m, y)

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
