
import numpy as np
from scipy.optimize import minimize_scalar

from chem_analysis.processing.processing_method import Baseline


class Subtract(Baseline):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 multiplier: float = 1,
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing, save_result)
        self.y_sub = y
        self.x_sub = x
        self.multiplier = multiplier

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if len(self.y_sub) == len(y):
            return self.multiplier * self.y_sub

        raise NotImplementedError()  # TODO: x-interpolation


class SubtractOptimize(Baseline):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 bounds: tuple[float, float] = (-2, 2),
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing, save_result)
        self.y_sub = y
        self.x_sub = x
        self.bounds = bounds

        # save
        self.multiplier = 1

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        self.multiplier = self._get_multiplier(x, y, self.x_sub, self.y_sub)
        return self.multiplier * self.y_sub

    def _get_multiplier(self, x: np.ndarray, y: np.ndarray, x_sub: np.ndarray, y_sub: np.ndarray) -> float:
        if len(self.y_sub) == len(y):
            def func(m) -> float:
                return float(np.sum(np.abs(y-m*y_sub)))

            result = minimize_scalar(func, bounds=self.bounds)
            if not result.success:
                raise ValueError(f"'{type(self).__name__}' has not converged.")
            return result.x

        raise NotImplementedError()  # TODO: x-interpolation
