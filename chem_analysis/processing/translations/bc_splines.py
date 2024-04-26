from typing import Iterable

import numpy as np
from scipy.optimize import minimize_scalar

from chem_analysis.processing.processing_method import ProcessingMethod
from chem_analysis.processing.weigths.weights import DataWeight, DataWeightChain


import numpy as np
from scipy.interpolate import UnivariateSpline

class Spline(BaselineCorrection):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 multiplier: float = 1,
                 ):
        super().__init__(None)
        self.y_sub = y
        self.x_sub = x
        self.multiplier = multiplier

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if len(self.y_sub) == len(y):
            # return self.multiplier * self.y_sub

        spline = UnivariateSpline(x, y, s=0)

        raise NotImplementedError()  # TODO: x-interpolation

    def get_baseline2D(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)

        for i in range(z.shape[0]):
            baseline[i, :] = self.get_baseline(x, z[i, :])

        return baseline
