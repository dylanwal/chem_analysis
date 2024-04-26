from typing import Iterable

import numpy as np

from chem_analysis.processing.processing_method import BaselineCorrection
from chem_analysis.processing.weigths.weights import DataWeight


class Polynomial(BaselineCorrection):
    def __init__(self,
                 degree: int = 1,
                 poly_weights: np.ndarray = None,
                 weights: DataWeight | Iterable[DataWeight] = None
                 ):
        super().__init__(weights)
        self.degree = degree
        self.poly_weights = poly_weights

    def run(self, x: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self.baseline = self.get_baseline(x, data)
        self.x = x
        self.data = data - self.baseline
        return x, self.data

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def get_baseline(self, x: np.ndarray, y: np.ndarray, poly_weights: np.ndarray = None) -> np.ndarray:
        if poly_weights is None:
            if self.poly_weights is None:
                poly_weights = np.ones_like(y)
            else:
                poly_weights = self.poly_weights

        if self.weights is not None:
            mask = self.weights.get_mask(x, y)
            x_ = x[mask]
            y_ = y[mask]
            w_ = poly_weights[mask]
        else:
            x_ = x
            y_ = y
            w_ = poly_weights

        params = np.polyfit(x_, y_, self.degree, w=w_)
        func_baseline = np.poly1d(params)
        return func_baseline(x)

    def get_baseline2D(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
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
