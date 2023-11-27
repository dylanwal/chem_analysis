import abc

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.utils.general_math import get_slice


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
    def __init__(self, y_subtract: np.ndarray, x_subtract: np.ndarray = None):
        self.y_subtract = y_subtract
        self.x_subtract = x_subtract  #TODO: add

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, y - self.y_subtract


class AlignMax(Translations):
    def __init__(self, range_: tuple[int, int] = None, range_index: slice = None):
        self.range_ = range_
        self.range_index = range_index
        self.max_ = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.range_ is not None:
            slice_ = get_slice(x, self.range_[0], self.range_[1])
            y_ = y[slice_]
        elif self.range_index is not None:
            y_ = y[self.range_index]
        else:
            y_ = y

        max_indices = np.argmax(y_)
        roll_amount = max_indices[0]

        x = np.roll(x, roll_amount)
        y = np.roll(y, roll_amount)
        return x, y
