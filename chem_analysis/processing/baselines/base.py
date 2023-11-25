import abc
from typing import Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.weigths.weights import DataWeights


class BaselineCorrection(ProcessingMethod, abc.ABC):
    def __init__(self, weights: DataWeights | Iterable[DataWeights] = None):
        self.weights = weights
        self._x = None
        self._y = None

    def get_weights(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.weights is None:
            return np.ones_like(y)
        if isinstance(self.weights, DataWeights):
            return self.weights.get_inverted_weights(x, y)

        weights = np.ones_like(y)
        for weight in self.weights:
            weights *= weight.get_inverted_weights(x, y*weights)
        return weights

    def apply_weights(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x, y = np.copy(x), np.copy(y)
        if self.weights is None:
            return x, y
        if isinstance(self.weights, DataWeights):
            return self.weights.apply_as_mask(x, y)

        for mask in self.weights:
            x, y = mask.apply_as_mask(x, y)

        return x, y

    @property
    def y(self) -> np.ndarray | None:
        return self._y

    @property
    def x(self) -> np.ndarray | None:
        return self._x

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class Polynomial(BaselineCorrection):
    def __init__(self, degree: int = 1, weights: DataWeights | Iterable[DataWeights] = None):
        super().__init__(weights)
        self.degree = degree

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        weights = self.get_weights(x, y)
        params = np.polyfit(x, y, self.degree, w=weights)
        func_baseline = np.poly1d(params)
        y_baseline = func_baseline(x)
        y = y - y_baseline

        self._y = y_baseline
        self._x = x
        return x, y





# Bernstein polynomial (order = 3)
# Splines
# multipoint
