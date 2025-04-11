import abc

import numpy as np

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


class CriteriaBase(MixinSubClassList, abc.ABC):
    @abc.abstractmethod
    def __call__(self, x: np.ndarray, y: np.ndarray, peak: np.ndarray, bound: np.ndarray) -> bool:
        ...


class CriteriaMaxPeaks(CriteriaBase):
    """
    Stop Criteria: Max Peaks

    Stop the algorithm after 'num_eval' iterations.
    * an iteration may be 1 function evaluation; but it may not be.

    """
    def __init__(self, max_peaks: int):
        self.max_peaks = max_peaks
        self._num_peaks = 0

    def __repr__(self):
        return f"{type(self).__name__} | max_peaks: {self.max_peaks}"

    def __call__(self, x: np.ndarray, y: np.ndarray, peak: np.ndarray, bound: np.ndarray) -> bool:
        self._num_peaks += 1
        if self._num_peaks >= self.max_peaks:
            return True
        return False


class CriteriaMinPeakHeight(CriteriaBase):
    """
    Peak height drops below threshold

    """
    def __init__(self, threshold: float | int, relative: bool = False):
        self.threshold = threshold
        self.relative = relative

    def __repr__(self):
        return f"{type(self).__name__} | threshold: {self.threshold}"

    def __call__(self, x: np.ndarray, y: np.ndarray, peak: np.ndarray, bound: np.ndarray) -> bool:
        threshold = np.max(y) * self.threshold if self.relative else self.threshold
        if threshold >= y[peak]:
            return True
        return False
