
import numpy as np
from scipy.optimize import minimize_scalar

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.base_obj.signal_ import Signal


class Subtract(Baseline):
    def __init__(self,
                 value: int | float | str = "min",
                 multiplier: int | float = 1,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """
        Subtract value from y.

        Parameters
        ----------
        value: int | float | str
            "min": subtract min value from y.
            "max": subtract max value from y.
            "mean": subtract average value from y.
            if int or float that value will be used.
        multiplier: int | float
            multiplier applied to value
        temporal_processing: int
            for 2d and 3d signals; they will be treated as a timeseries and each row will be processed as 1 dimensional
            lower signal.
        save_result:
            True: saves resulting baseline in 'Baseline.baseline' (typically done for plotting)
            False: result baseline not saved. (typically done to minimize memory usage)
        """
        super().__init__(temporal_processing, save_result)
        self.value = None
        if isinstance(value, str):
            if not (value == "min" or value == "max"):
                raise ValueError(f"'value' must be either 'min', 'max', or numerical value. \n\tGiven 'value': {value}")

        self._value = value
        self.multiplier = multiplier

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self._value == "min":
            self.value = np.min(y)
        if self._value == "max":
            self.value = np.max(y)
        if self._value == "mean":
            self.value = np.mean(y)

        return np.ones_like(y) * self.value * self.multiplier


class SubtractSignal(Baseline):
    def __init__(self,
                 signal: Signal,
                 multiplier: float = 1,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """
        Subtract signal from y.

        Parameters
        ----------
        signal: Signal
            signal to be used for subtraction
        multiplier: int | float
            multiplier applied to value
        temporal_processing: int
            for 2d and 3d signals; they will be treated as a timeseries and each row will be processed as 1 dimensional
            lower signal.
        save_result:
            True: saves resulting baseline in 'Baseline.baseline' (typically done for plotting)
            False: result baseline not saved. (typically done to minimize memory usage)
        """
        super().__init__(temporal_processing, save_result)
        self.signal = signal
        self.multiplier = multiplier

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if len(self.signal.y) != len(y):
            y_sub = np.interp(x, self.signal.x, self.signal.y)
        else:
            y_sub = self.signal.y

        return self.multiplier * y_sub


class SubtractSignalOptimize(Baseline):
    def __init__(self,
                 signal: Signal,
                 bounds: tuple[float, float] = (-2, 2),
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        """
        Subtract signal from y with optimization. The optimization will find max multiplayer as to not
        create negative values.

        Parameters
        ----------
        signal: Signal
            signal to be used for subtraction
        bounds: tuple[float, float]
            range that the multiplier will explore
        temporal_processing: int
            for 2d and 3d signals; they will be treated as a timeseries and each row will be processed as 1 dimensional
            lower signal.
        save_result:
            True: saves resulting baseline in 'Baseline.baseline' (typically done for plotting)
            False: result baseline not saved. (typically done to minimize memory usage)
        """
        super().__init__(temporal_processing, save_result)
        self.signal = signal
        self.bounds = bounds

        # save
        self.multiplier = 1

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if len(self.signal.y) != len(y):
            y_sub = np.interp(x, self.signal.x, self.signal.y)
        else:
            y_sub = self.signal.y

        self.multiplier = self._get_multiplier(y, y_sub)
        return self.multiplier * y_sub

    def _get_multiplier(self, y: np.ndarray, y_sub: np.ndarray) -> float:
        def func(m) -> float:
            offset = y-m*y_sub
            mask = offset > 0
            return float(np.sum(offset[mask]) - np.sum(offset[np.logical_not(mask)]*10))        # float(np.sum(np.abs(y-m*y_sub)))

        result = minimize_scalar(func, bounds=self.bounds)
        if not result.success:
            raise ValueError(f"'{type(self).__name__}' has not converged.")
        return result.x
