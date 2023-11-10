import abc

import numpy as np

from chem_analysis.processing.base import ProcessingMethod


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


class Subtract(Translations):
    def __init__(self, y_subtract: np.ndarray):
        self.y_subtract = y_subtract

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, y - self.y_subtract


class AlignMax(Translations):
    def __init__(self, range_: tuple[int, int] = None, range_index: tuple[int, int] = None):
        self.range_ = range_
        self.range_index = range_index
        self.max_ = None

    def align_max(self, x: np.ndarray, y: np.ndarray):
        max_indices = np.argmax(data, axis=1)
        roll_amount = -max_indices + max_indices[0]
        for i in range(1, data.shape[0]):
            data[i] = np.roll(data[i], roll_amount[i])
        return data