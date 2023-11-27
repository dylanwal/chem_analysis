import abc
from typing import Iterable, Callable, Sequence

import numpy as np

from chem_analysis.utils.general_math import get_slice
from chem_analysis.utils.code_for_subclassing import MixinSubClassList
import chem_analysis.processing.weigths.penalty_functions as penalty_functions


class DataWeights(MixinSubClassList, abc.ABC):
    def __init__(self, threshold: float = 0.5, normalized: bool = True):
        self._weights = None
        self.threshold = threshold
        self.normalized = normalized

    @property
    def weights(self) -> np.ndarray | None:
        return self._weights

    @abc.abstractmethod
    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        ...

    def get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self._weights is None:
            self._weights = self._get_weights(x, y)
            if np.all(self._weights == 0):
                raise ValueError(f"All weights are zero after applying {type(self).__name__}")

        return self._weights

    def get_normalized_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        weights = self.get_weights(x, y)
        return weights / np.max(weights)

    def get_inverted_weights(self, x: np.ndarray, y: np.ndarray, replace_zero: float = None) -> np.ndarray:
        weights = self.get_weights(x, y)

        # avoid divide by zero
        if replace_zero is None:
            replace_zero = np.min(weights[weights > 0]) * 0.9
            # 0.9 is just to make it just a bit smaller than the smallest value

        mask = weights == 0
        weights[mask] = replace_zero

        return 1 / weights

    def apply_as_mask(self, x: np.ndarray, y: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray]:
        """ Only use first row of mask because not returning multiple x. """
        if self.normalized:
            weights = self.get_normalized_weights(x, y)
        else:
            weights = self.get_weights(x, y)

        indexes = weights >= self.threshold
        if len(y.shape) == 1:
            return x[indexes], y[indexes]
        return x[indexes[0, :]], y[:, indexes[0, :]]


class Slices(DataWeights):
    def __init__(self,
                 slices: slice | Iterable[slice],
                 threshold: float = 0.5,
                 normalized: bool = True,
                 invert: bool = False
                 ):
        super().__init__(threshold, normalized)
        self.slices = slices
        self.invert = invert

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if not isinstance(self.slices, list):
            slices = [self.slices]
        else:
            slices = self.slices

        weights = np.ones_like(y)

        for slice_ in slices:
            if len(y.shape) == 1:
                weights[slice_] = 0
            else:
                weights[:, slice_] = 0

        if self.invert:
            return np.logical_not(weights)
        return weights


class Spans(DataWeights):
    def __init__(self,
                 x_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2
                 threshold: float = 0.5,
                 normalized: bool = True,
                 invert: bool = False
                 ):
        super().__init__(threshold, normalized)
        self.x_spans = x_spans
        self.invert = invert

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if not isinstance(self.x_spans, list):
            x_spans = [self.x_spans]
        else:
            x_spans = self.x_spans

        weights = np.ones_like(y)

        for x_span in x_spans:
            slice_ = get_slice(x, x_span[0], x_span[1])
            if len(y.shape) == 1:
                weights[slice_] = 0
            else:
                weights[:, slice_] = 0

        if self.invert:
            return np.logical_not(weights)
        return weights


class MultiPoint(DataWeights):
    def __init__(self, indexes: Iterable[int], threshold: float = 0.5, normalized: bool = True):
        """

        Parameters
        ----------
        indexes:
            index where values will be kept
        """
        super().__init__(threshold, normalized)
        self.indexes = indexes

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        weights = np.zeros_like(y)
        weights[self.indexes] = 1
        return weights


class DistanceMedian(DataWeights):
    def __init__(self,
                 penalty_function: Callable = penalty_functions.penalty_function_linear,
                 threshold: float = 0.5,
                 normalized: bool = True
                 ):
        super().__init__(threshold, normalized)
        self.penalty_function = penalty_function

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return self.penalty_function(y-np.median(y))


def get_distance_remove(y: np.ndarray, amount: float = 0.7, speed: float = 0.1, max_iter: int = 1000):
    stop_len = len(y) * amount
    indexes = np.ones_like(y, dtype=bool)
    for i in range(max_iter):
        dist_from_median = np.abs(y - np.median(y[indexes]))
        median_deviation = np.median(dist_from_median[indexes])
        if median_deviation == 0:
            median_deviation = np.mean(dist_from_median[indexes])
        scale_distances_from_median = dist_from_median / median_deviation
        cut_off_distance = np.max(scale_distances_from_median[indexes]) * (1 - 0.5 * np.exp(-i / speed))
        keep_points = scale_distances_from_median < cut_off_distance
        if np.sum(keep_points) == 0:
            break
        indexes = np.bitwise_and(indexes, keep_points)

        if np.sum(indexes) < stop_len:
            break

    return indexes


class AdaptiveDistanceMedian(DataWeights):
    def __init__(self,
                 penalty_function: Callable = penalty_functions.penalty_function_linear,
                 amount: float = 0.5,
                 speed: float = 10,
                 max_iter: int = 1000,
                 threshold: float = 0.5,
                 normalized: bool = True
                 ):
        super().__init__(threshold, normalized)
        self.amount = amount
        self.speed = speed
        self.max_iter = max_iter
        self.indexes = None
        self.penalty_function = penalty_function

    def _get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        self.indexes = get_distance_remove(y, self.amount, self.speed, self.max_iter)
        return self.penalty_function(y-np.median(y[self.indexes]))
