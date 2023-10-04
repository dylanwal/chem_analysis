import abc

import numpy as np

from chem_analysis.analysis.algorithms.processing.base import ProcessingMethod


class Translations(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class Horizontal(Translations):
    def __init__(self, shift: float = 1, wrap: bool = False):
        self.shift = shift  # Hz

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        y = y[self.shiftPoints:]
        nmrData.fidTimeForLB = nmrData.fidTime[self.shiftPoints:]
        return x, y

