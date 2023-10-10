import abc

import numpy as np
from scipy.fft import fft
from scipy.fftpack import fftshift

from chem_analysis.algorithms.processing.base import ProcessingMethod


class FourierTransform(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class FastFourierTransform(FourierTransform):
    def __init__(self, degree: int = 1):
        self.degree = degree
        self._y_baseline = None

    @property
    def y_baseline(self) -> np.ndarray:
        return self._y_baseline

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, fftshift(fft(y))
