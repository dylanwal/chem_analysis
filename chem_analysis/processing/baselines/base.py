import abc
from typing import Iterable

import numpy as np

from chem_analysis.processing.base import ProcessingMethod
from chem_analysis.processing.weigths.weights import DataWeight, DataWeightChain, DataWeightsBase


class BaselineCorrection(ProcessingMethod, abc.ABC):
    def __init__(self, weights: DataWeight | Iterable[DataWeight] | DataWeightChain = None):
        if weights is not None and not isinstance(weights, DataWeight) and not isinstance(weights, DataWeightChain):
            weights = DataWeightChain(weights)
        self.weights: DataWeightsBase = weights
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
    def get_baseline(self, x: np.ndarray, y: np.ndarray, weights: np.ndarray = None) -> np.ndarray:
        ...

    def get_baseline_array(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)
        if self.weights:
            x, y, z = self.weights.apply_as_mask_array_index(x, y, z)

        for i in range(z.shape[0]):
            baseline[i, :] = self.get_baseline(x, z[i, :])

        return baseline


class Polynomial(BaselineCorrection):
    def __init__(self, degree: int = 1, weights: DataWeight | Iterable[DataWeight] = None):
        super().__init__(weights)
        self.degree = degree

    def get_baseline(self, x: np.ndarray, y: np.ndarray, weights: np.ndarray = None) -> np.ndarray:
        x_, y_ = self.weights.apply_as_mask(x, y)
        if weights is None:
            weights = np.ones_like(y_)

        params = np.polyfit(x_, y_, self.degree, w=weights)
        func_baseline = np.poly1d(params)
        return func_baseline(x)


# Bernstein polynomial (order = 3)
# Splines
# multipoint
