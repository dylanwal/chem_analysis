from typing import Callable

import numpy as np
import scipy.optimize


class Calibration:
    def __init__(self,
                 calibration_function: Callable,
                 lower_bound_y: int | float = None,
                 upper_bound_y: int | float = None,
                 name: str = None
                 ):
        self.name = name
        self.calibration_function = calibration_function
        self.lower_bound_y = lower_bound_y
        self.upper_bound_y = upper_bound_y

        self._upper_bound_x = None
        self._lower_bound_x = None

    @property
    def upper_bound_location(self) -> None | float:
        if self.upper_bound_y is None:
            return None
        if self._upper_bound_x is None:
            def _ub(x: np.ndarray) -> np.ndarray:
                return self.get_y(x) - self.upper_bound_y

            result_ub = scipy.optimize.root_scalar(_ub, x0=0.1, x1=0.2)
            self._upper_bound_location = result_ub.root

        return self._upper_bound_location

    @property
    def lower_bound_location(self) -> None | float:
        if self.lower_bound_y is None:
            return None

        if self._lower_bound_x is None:
            def _lb(x: np.ndarray) -> np.ndarray:
                return self.get_y(x) - self.lower_bound_y

            if self.upper_bound_location:
                x0 = self.upper_bound_location
            else:
                x0 = 0.1

            result_lb = scipy.optimize.root_scalar(_lb, x0=x0, x1=x0 + 0.1)
            self._lower_bound_location = result_lb.root

        return self._lower_bound_location

    def check_bounds(self, x: int | float | np.ndarray, error: bool = False):
        # TODO

    def get_y(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if self.lower_bound_x is not None:
            if

        return self.calibration_function(x)
