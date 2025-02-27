import logging

import numpy as np
from scipy.ndimage import grey_opening, grey_erosion, grey_dilation

from chem_analysis.processing.processing_method import Baseline

logger = logging.getLogger(__name__)


def norm2(old: np.ndarray, new: np.ndarray) -> float:
    return np.linalg.norm(new - old) / np.linalg.norm(old)


def average_opening(y: np.ndarray, window_size: int = 10, opening: bool = True) -> np.ndarray:
    if not opening:
        y = average_opening(y, window_size)
    return (grey_dilation(y, window_size) + grey_erosion(y, window_size)) / 2


def morphological_average(y: np.ndarray, window_size: int = 10) -> np.ndarray:
    opening = grey_opening(y, window_size)
    return np.minimum(opening, average_opening(opening, window_size))


def iterative_morphological_baseline(y: np.ndarray, window_size: int = 10, error: float = 1e-4,
                                     max_iter: int = 100) -> np.ndarray:
    baseline = y
    for i in range(max_iter):
        new_baseline = morphological_average(y, window_size)
        if norm2(baseline, new_baseline) < error:
            return baseline
        baseline = new_baseline

    logger.warning('Maximum iterations exceeded prior to finding a solution within error tolerance.')
    return baseline


def estimate_window(
        y: np.ndarray,
        min_size: int = 3,
        max_size: int | None = None,
        steps: int = 30,
        error: float | int = 1e-5,
        n: int = 3
) -> int:
    """
    Iteratively change window size till 'n' consecutive widows are equal.
    https://doi.org/10.1177/0003702817752371

    Parameters
    ----------
    y
    min_size
    max_size
    steps
    error
    n

    Returns
    -------

    """
    if max_size is None:
        max_size = (y.size - 1) // 2

    old = grey_opening(y, min_size)
    windows = np.linspace(min_size, max_size, steps, dtype=int)

    best_window, count = min_size, 0
    for i in range(steps):  # TODO: replace grid search with optimization
        new = grey_opening(y, windows[i])
        if norm2(old, new) < error:
            if count > n:
                return int(windows[i - n])
            else:
                count += 1

        else:
            count = 0

    return int(windows[int(len(windows) / 3)])


class MorphologicalAverage(Baseline):
    def __init__(self,
                 window_size: int = None,
                 invert: bool = False,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """

        https://doi.org/10.1177/0003702817752371

        Parameters
        ----------
        window_size: int
            size of window

        Returns
        -------
        mask:
            1 where the baseline is
            0 where peaks are
    """
        super().__init__(temporal_processing, save_result)
        self.window = window_size

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window is None:
            self.window = estimate_window(y)
        return morphological_average(y, self.window)


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

        Returns
        -------
        mask:
            1 where the baseline is
            0 where peaks are
    """
        super().__init__(temporal_processing, save_result)
        self.window = None
        self.max_iter = max_iter
        self.error = error

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.window is None:
            self.window = estimate_window(y)
        return iterative_morphological_baseline(y, self.window, self.error, self.max_iter)
