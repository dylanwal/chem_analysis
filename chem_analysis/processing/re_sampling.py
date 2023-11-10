import abc

import numpy as np

from chem_analysis.processing.base import ProcessingMethod


class ReSampling(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class EveryN(ReSampling):
    def __init__(self, n: int = 2):
        self.n = n
        self._y_baseline = None

    @property
    def y_baseline(self) -> np.ndarray:
        return self._y_baseline

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x[::self.n], y[::self.n]


# class FitResample(ReSampling):
#     def __init__(self, n: int = 2):
#         self.n = n
#         self._y_baseline = None
#
#     @property
#     def y_baseline(self) -> np.ndarray:
#         return self._y_baseline
#
#     def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
#         return x[::self.n], y[::self.n]