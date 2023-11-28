import abc
from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.utils.general_math import get_slice


class Remove(ProcessingMethod, abc.ABC):
    def __init__(self):
        self.index = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError("Signals can't remove signals")
        return None  # noqa

    def run_array(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.index = self.get_index(x, y, z)
        return x, y[self.index], z[self.index]

    @abc.abstractmethod
    def get_index(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        ...


class EveryN(Remove):
    def __init__(self, n: int = 2, start_index: int = 0):
        super().__init__()
        self.n = n
        self.start_index = start_index

    def get_index(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        index = np.ones_like(y, dtype=np.bool)
        index[self.start_index::self.n, :] = 0
        return index


class CutOffValue(Remove):
    def __init__(self,
                 x_span: float | Sequence[float],
                 cut_off_value: float | int,
                 invert: bool = False
                 ):
        super().__init__()
        self.x_span = x_span
        self.cut_off_value = cut_off_value
        self.invert = invert

    def get_index(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        if isinstance(self.x_span, Sequence):
            slice_ = get_slice(x, self.x_span[0], self.x_span[1])
            indexes = np.any(z[:, slice_] > self.cut_off_value)
        else:
            index = np.argmin(np.abs(x - self.x_span))
            indexes = z[:, index] < self.cut_off_value

        if self.invert:
            return np.logical_not(indexes)
        return indexes
