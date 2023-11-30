import abc
from typing import Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.weigths.weights import DataWeight, DataWeightChain


class BaselineCorrection(ProcessingMethod, abc.ABC):
    def __init__(self, weights: DataWeight | Iterable[DataWeight] = None):
        if weights is not None and isinstance(weights, Iterable):
            weights = DataWeightChain(weights)
        self.weights: DataWeight = weights
        self._x = None
        self._y = None

    @property
    def y(self) -> np.ndarray | None:
        return self._y

    @property
    def x(self) -> np.ndarray | None:
        return self._x

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self._y = self.get_baseline(x, y)
        self._x = x
        return x, y - self._y

    def run_array(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._y = self.get_baseline_array(x, y, z)
        self._x = x
        return x, y, z - self._y

    @abc.abstractmethod
    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        ...

    def get_baseline_array(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)

        for i in range(z.shape[0]):
            baseline[i, :] = self.get_baseline(x, z[i, :])

        return baseline


class Polynomial(BaselineCorrection):
    def __init__(self,
                 degree: int = 1,
                 poly_weights: np.ndarray = None,
                 weights: DataWeight | Iterable[DataWeight] = None
                 ):
        super().__init__(weights)
        self.degree = degree
        self.poly_weights = poly_weights

    def get_baseline(self, x: np.ndarray, y: np.ndarray, poly_weights: np.ndarray = None) -> np.ndarray:
        if poly_weights is None:
            if self.poly_weights is None:
                poly_weights = np.ones_like(y)
            else:
                poly_weights = self.poly_weights

        mask = self.weights.get_mask(x, y)

        params = np.polyfit(x[mask], y[mask], self.degree, w=poly_weights[mask])
        func_baseline = np.poly1d(params)
        return func_baseline(x)

    def get_baseline_array(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)

        if self.poly_weights is None:
            self.poly_weights = np.ones_like(z[0, :])

        if self.poly_weights.shape == z.shape:
            for i in range(z.shape[0]):
                baseline[i, :] = self.get_baseline(x, z[i, :], self.poly_weights[i, :])
        elif self.poly_weights.size == z.shape[1]:
            for i in range(z.shape[0]):
                baseline[i, :] = self.get_baseline(x, z[i, :], self.poly_weights)
        else:
            raise ValueError(f"{type(self).__name__}.poly_weights is wrong shape."
                             f"\n\texpected: {z.shape} or {z.shape[1]}"
                             f"\n\tgiven: {self.poly_weights.shape}")

        return baseline


# Bernstein polynomial (order = 3)
# Splines
# multipoint
