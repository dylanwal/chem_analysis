import abc
from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.weigths.weights import Slices, Spans


class ReSampling(ProcessingMethod, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class EveryN(ReSampling):
    def __init__(self, n: int = 2):
        self.n = n

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x[::self.n], y[::self.n]


class CutSlices(ReSampling):
    def __init__(self,
                 slices: slice | Iterable[slice],
                 invert: bool = True
                 ):
        self.slices = slices
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        span = Slices(self.slices)
        mask = span.apply_as_mask(x, y)

        if len(y.shape) == 1:
            return x[mask], y[mask]
        return x[mask], y[:, mask]


class CutSpans(ReSampling):
    def __init__(self,
                 x_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2
                 invert: bool = True
                 ):
        self.x_spans = x_spans
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        span = Spans(self.x_spans)
        mask = span.apply_as_mask(x, y)

        if len(y.shape) == 1:
            return x[mask], y[mask]
        return x[mask], y[:, mask]


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