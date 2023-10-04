import abc

import numpy as np

from chem_analysis.analysis.algorithms.processing.base import ProcessingMethod


class Smoothing(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class LineBroadening(Smoothing):
    def __init__(self, degree: float = 1):
        self.degree = degree  # Hz

    @property
    def y_baseline(self) -> np.ndarray:
        return self._y_baseline

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        length = len(nmrData.allFid[-1][0])
        sp.multiply(nmrData.allFid[-1][:], sp.exp(-nmrData.fidTimeForLB[:length] * self.degree * np.pi))
        return x, y

