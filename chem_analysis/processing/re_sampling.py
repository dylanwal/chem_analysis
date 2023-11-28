import abc
from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.weigths.weights import Slices, Spans


class ReSampling(ProcessingMethod, abc.ABC):
    ...


class EveryN(ReSampling):
    def __init__(self, n: int = 2, start_index: int = 0):
        self.n = n
        self.start_index = start_index

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x[::self.n], y[::self.n]

    def run_array(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x[self.start_index::self.n], y, z[:, self.start_index::self.n]


class CutSlices(ReSampling):
    def __init__(self,
                 slices: slice | Iterable[slice],
                 invert: bool = False
                 ):
        self.slices = slices
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        slice_ = Slices(self.slices, invert=self.invert)
        return slice_.apply_as_mask(x, y)

    def run_array(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        slice_ = Slices(self.slices, invert=self.invert)
        x, z = slice_.apply_as_mask(x, z)
        return x, y, z


class CutSpans(ReSampling):
    def __init__(self,
                 x_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2
                 invert: bool = False
                 ):
        self.x_spans = x_spans
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        span = Spans(self.x_spans, invert=self.invert)
        return span.apply_as_mask(x, y)

    def run_array(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        span = Spans(self.x_spans, invert=self.invert)
        x, z = span.apply_as_mask(x, z)
        return x, y, z

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