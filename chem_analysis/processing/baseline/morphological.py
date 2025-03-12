import logging

import numpy as np
from scipy.ndimage import grey_opening, grey_erosion, grey_dilation, grey_closing, minimum_filter
from scipy.signal import convolve

from chem_analysis.processing.processing_method import Baseline
import chem_analysis.utils.math as math_utils

logger = logging.getLogger(__name__)


def average_opening(y: np.ndarray, window_size: int = 10, opening: bool = True) -> np.ndarray:
    if not opening:
        y = average_opening(y, window_size)
    return (grey_dilation(y, window_size) + grey_erosion(y, window_size)) / 2


def morphological_average(y: np.ndarray, window_size: int = 10) -> np.ndarray:
    if len(y)/window_size < 3:
        logger.warning('"window_size" may be too large causing issues. ')
    opening = grey_opening(y, window_size)
    return np.minimum(opening, average_opening(opening, window_size))


def average_open_close(y: np.ndarray, window_size: int = 10) -> np.ndarray:
    return (grey_closing(y, window_size) + grey_opening(y, window_size))/2


def iterative_morphological_baseline(
        y: np.ndarray,
        window_size: int = 10,
        error: float = 1e-4,
        max_iter: int = 100
) -> np.ndarray:
    baseline = y
    for i in range(max_iter):
        new_baseline = morphological_average(y, window_size)
        if math_utils.norm2(baseline, new_baseline) < error:
            return baseline
        baseline = new_baseline

    logger.warning('Maximum iterations exceeded prior to finding a solution within error tolerance.')
    return baseline


def iterative_morphological_mollify(y: np.ndarray, window_size: int = 10, num_iter: int = 100, error: float = 1e-6):
    """
    some testing suggest there may be an issue; the error is not good stopping condtion and num_iter changes
    result quite a bit
    so not adding to list at this point in time
    https://doi.org/10.1177/0003702818811688
    """
    mollify_base = mollify_kernel(window_size)
    conv_edge_base = convolve(np.ones_like(y), mollify_base, mode='same')
    window_size = int(2*window_size+1)
    baseline = np.copy(y)
    for k in range(num_iter):
        b = np.minimum(y, average_open_close(baseline, window_size))
        new_baseline = convolve(b, mollify_base, mode='same') / conv_edge_base
        if k >= 2 and math_utils.norm2(baseline, new_baseline) <= error:
            break
        baseline = new_baseline

    return baseline


def mollify_kernel(window_size: int) -> np.ndarray:
    """
    Applies the Gaussian mollifier to a 1D signal.
    """
    x = (np.arange(0, 2 * window_size + 1) - window_size) / window_size
    kernel = np.zeros_like(x)
    kernel[1:-1] = np.exp(-1 / (1 - (x[1:-1]) ** 2))
    return kernel / kernel.sum()


def morphological_mollifer(
        y: np.ndarray,
        smooth_width: int = 6,
        window_size: int = 50,
        max_iter: int = 25,
        error: float = 1e-1
) -> np.ndarray:
    # # Define mollifier kernels
    mollify_spec = mollify_kernel(smooth_width)
    conv_edge_spec = convolve(np.ones_like(y), mollify_spec, mode='same')
    mollify_base = mollify_kernel(window_size)
    conv_edge_base = convolve(np.ones_like(y), mollify_base, mode='same')

    baseline_smooth = np.zeros((len(y), max_iter))
    for k in range(max_iter):
        ysmooth = convolve(y, mollify_spec, mode='same') / conv_edge_spec
        baseline = minimum_filter(ysmooth, size=window_size, mode='nearest')

        baseline_smooth[:, k] = convolve(baseline, mollify_base, mode='same') / conv_edge_base
        if k >= 2 and math_utils.norm2(baseline_smooth[:, k - 1], baseline_smooth[:, k]) <= error:
            baseline_smooth = baseline_smooth[:, :k]
            print(k)
            break
        y -= baseline_smooth[:, k]

    logger.warning('Maximum iterations exceeded prior to finding a solution within error tolerance.')
    return np.sum(baseline_smooth, axis=1)


def estimate_window(
        y: np.ndarray,
        min_size: int | None = None,
        max_size: int | None = None,
        error: float | int = 1e-5,
        n: int = 3
) -> int:
    """
    Iteratively change window size till 'n' consecutive widows are equal.
    https://doi.org/10.1177/0003702817752371
    DOI: 10.1366/000370210791414281

    The optimal size of the structuring element is thought to be at least equal to or slightly larger
    than the widest peak feature in the original signal.

    Parameters
    ----------
    y
    min_size
    max_size
    error
    n

    Returns
    -------

    """
    if min_size is None:
        min_size = 3
    if max_size is None:
        max_size = y.size // 2

    old = grey_opening(y, min_size)
    window, count = min_size, 0
    for i in range((max_size-min_size)//2):
        window += 2
        new = grey_opening(y, window)
        if math_utils.norm2(old, new) < error:
            if count >= n:
                return int((window - 6)/2) + 1
            else:
                count += 1

        else:
            count = 0
        old = new

    return min([3, int((max_size-min_size) / 20)])  # pick something at a 1/10 if the search fails / peaks assumed to be small in signal


class MorphologicalAverage(Baseline):
    def __init__(self,
                 window_size: int = None,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        https://doi.org/10.1177/0003702817752371

        Parameters
        ----------
        window_size: int
            size of window

    """
        super().__init__(temporal_processing, save_result)
        self.window_size = window_size

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window_size is None:
            self.window_size = estimate_window(y)
        return morphological_average(y, self.window_size)


class MorphologicalAutoWindow(Baseline):
    def __init__(self,
                 max_iter: int = 100,
                 error: float = 1e-4,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        https://doi.org/10.1177/0003702817752371

        Parameters
        ----------
        max_iter:
            max_iterations to refine baseline

    """
        super().__init__(temporal_processing, save_result)
        self.window_size = None
        self.max_iter = max_iter
        self.error = error

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window_size is None:
            self.window_size = estimate_window(y)
        return iterative_morphological_baseline(y, self.window_size, self.error, self.max_iter)


class MorphologicalMollifier(Baseline):
    def __init__(self,
                 smooth_width: int = None,
                 window_size: int = None,
                 error: int | float = 1e-1,
                 max_iter: int = 25,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

         https://doi.org/10.1002/jrs.5010

        Parameters
        ----------
        smooth_width:
            higher number does more smoothing
        window_size:
            size of window
        max_iter:
            max_iterations to refine baseline

    """
        super().__init__(temporal_processing, save_result)
        self.window_size = window_size
        self.smooth_width = smooth_width
        self.max_iter = max_iter
        self.error = error

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window_size is None:
            self.window_size = estimate_window(y)
        if self.smooth_width is None:
            self.smooth_width = max([1, int(self.window_size/4)])  # guess; probably better methods
        return iterative_morphological_baseline(y, self.window_size, self.error, self.max_iter)


class MorphologicalTopHat(Baseline):
    def __init__(self,
                 window_size: int = None,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        Parameters
        ----------

        window_size:
            size of window

    """
        super().__init__(temporal_processing, save_result)
        self.window_size = window_size

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window_size is None:
            self.window_size = estimate_window(y)

        return grey_opening(y, self.window_size)


def distance_to_circle(y: float, x_y: np.ndarray, center_index: int, radius: float) -> float:
    """Calculate the minimum distance from points to an adjusted circle center."""
    center_circle = np.array([x_y[center_index, 0], y])
    distances = np.linalg.norm(x_y - center_circle, axis=1)
    return np.min(distances) - radius


from scipy.optimize import toms748


def rolling_ball_baseline(x: np.ndarray, y: np.ndarray, radius: float = 3) -> np.ndarray:
    n = np.argmin(np.abs((x - x[0]) - 3))
    n = 2 * (n // 2) + 1  # make it always odd

    baseline = np.zeros_like(y)
    for i in range(len(y)):
        x_ = x[max(0, i - n):min(len(y), i + n)]
        y_ = y[max(0, i - n):min(len(y), i + n)]
        x_y = np.column_stack((x_, y_))
        center_index = (len(x_) // 2)
        result = toms748(distance_to_circle, np.min(y_) - 1.1 * radius, y_[center_index] - 0.9 * radius,
                         args=(x_y, center_index, radius))
        baseline[i] = result + radius

    return baseline


class MorphologicalBallRolling(Baseline):
    def __init__(self,
                 radius: int | float,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        Parameters
        ----------

        radius:
            size of ball to be rolled

    """
        super().__init__(temporal_processing, save_result)
        self.radius = radius

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return rolling_ball_baseline(x, y, self.radius)


def adaptive_ball_rolling(
        x: np.ndarray,
        y: np.ndarray,
        min_r: float | int | None = None,
        max_r: float = None,
        max_iter: int = 100,
        min_slope: float = 1e-4,
) -> np.ndarray:
    """
Adaptive Rolling Ball Baseline Correction with Bidirectional Smoothing.
"""
    if min_r is None:
        min_r = abs(x[0]-x[int(min([1, len(x)//500]))])
    if max_r is None:
        max_r = np.max(x) * 2

    r = np.linspace(min_r, max_r, max_iter, dtype=float)
    norms = np.zeros(max_iter)
    prev_baseline = np.zeros_like(y)
    for i in range(max_iter):
        baseline = rolling_ball_baseline(x, y, r[i])
        norms[i] = math_utils.norm2(prev_baseline, baseline)
        prev_baseline = baseline
        if i <= 3:
            continue
        slope_norms = math_utils.get_slope(norms[i-3:i])
        if abs(slope_norms) < min_slope:
            break

    else:
        print("did not converge")

    return prev_baseline


class MorphologicalAdaptiveBallRolling(Baseline):
    def __init__(self,
                min_r: float | int | None = None,
                max_r: float = None,
                max_iter: int = 100,
                min_slope: float = 1e-4,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        Parameters
        ----------

        radius:
            size of ball to be rolled

    """
        super().__init__(temporal_processing, save_result)
        self.min_r = min_r
        self.max_r = max_r
        self.max_iter = max_iter
        self.min_slope = min_slope

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return adaptive_ball_rolling(x, y, self.min_r, self.max_r, max_iter=self.max_iter, min_slope=self.min_slope)