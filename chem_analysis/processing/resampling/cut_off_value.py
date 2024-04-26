from typing import Sequence

import numpy as np

from chem_analysis.processing.processing_method import ReSampling
from chem_analysis.utils.math import get_slice


class CutOffValue(ReSampling):
    def __init__(self,
                 x_span: float | Sequence[float],  # TODO: generalize to n dimensions
                 cut_off_value: float | int,
                 invert: bool = False,
                 non_temporal_processing: bool = False
                 ):
        super().__init__(non_temporal_processing)
        self.x_span = x_span
        self.cut_off_value = cut_off_value
        self.invert = invert
        self.index = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError("Not valid method for Signals")

    def run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        index = self.get_index(x, z)
        self.index = index
        return x, y[index], z[index]

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def get_index(self, x: np.ndarray, z: np.ndarray) -> np.ndarray:
        if isinstance(self.x_span, Sequence):
            slice_ = get_slice(x, self.x_span[0], self.x_span[1])
            indexes = np.any(z[:, slice_] > self.cut_off_value)
        else:
            index = np.argmin(np.abs(x - self.x_span))
            indexes = z[:, index] < self.cut_off_value

        if self.invert:
            return np.logical_not(indexes)
        return indexes


