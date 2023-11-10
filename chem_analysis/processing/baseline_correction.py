import abc

import numpy as np

from chem_analysis.processing.base import ProcessingMethod


class BaselineCorrection(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class Polynomial(BaselineCorrection):
    def __init__(self, degree: int = 1):
        self.degree = degree
        self._y_baseline = None

    @property
    def y_baseline(self) -> np.ndarray:
        return self._y_baseline

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        params = np.polyfit(x, y, self.degree)
        func_baseline = np.poly1d(params)
        y_baseline = func_baseline(x)
        y = y - y_baseline

        self._y_baseline = y_baseline
        return x, y

# Whittaker smoother (filter=21 Hz, smoothfactor = 32768)
# Bernstein polynomial (order = 3)
# Ablative (points 5, passes 10)
# Splines
# multipoint