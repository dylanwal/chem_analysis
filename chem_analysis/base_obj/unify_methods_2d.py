import abc
from typing import Sequence

import numpy as np

import chem_analysis.utils.math as ca_math
from chem_analysis.base_obj.signal_2d import Signal2D


class UnifyMethod2D(abc.ABC):
    @abc.abstractmethod
    def run(self, data: Sequence[Signal2D]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...


class UnifyMethodStrict2D(UnifyMethod2D):
    def __init__(self, rtol: float | None = None, atol: float | None = None):
        self.rtol= rtol
        self.atol = atol

    def get_args(self) -> dict:
        if self.rtol is None and self.atol is None:
            self.rtol = 0.01
        dict_ = {}
        if self.rtol is not None:
            dict_['rtol'] = self.rtol
        if self.atol is not None:
            dict_['atol'] = self.atol
        return dict_

    def run(self, signals: Sequence[Signal2D]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x = signals[0].x
        y = signals[0].y
        z = np.empty((len(signals), len(y), len(x)), dtype=signals[0].z.dtype)
        args = self.get_args()
        for i, sig in enumerate(signals):
            if not np.all(np.isclose(sig.x, x, **args)):
                raise ValueError(f"Signal {i} has a different x-axis than first signal.")
            if not np.all(np.isclose(sig.y, y, **args)):
                raise ValueError(f"Signal {i} has a different y-axis than first signal.")
            z[i, :, :] = sig.z

        return x, y, z


class UnifyMethodMS2D(UnifyMethod2D):
    """
    ms values will be set as integers
    """
    def __init__(self, ms_max: int | None = None, ms_min: int | None = None):
        self.ms_max = ms_max
        self.ms_min = ms_min

    def get_x(self, signals: Sequence[Signal2D]) -> np.ndarray:
        if self.ms_max is None:
            max_ = np.max([sig.x[-1] for sig in signals])
        else:
            max_ = self.ms_max
        if self.ms_min is None:
            min_ = np.min([sig.x[0] for sig in signals])
        else:
            min_ = self.ms_min

        return np.arange(min_, max_ + 1, dtype=ca_math.min_uint_dtype(max_))

    def get_y(self, signals: Sequence[Signal2D]) -> np.ndarray:
        min_fid_length = min([len(sig.y_raw) for sig in signals])
        return signals[0].y_raw[:min_fid_length]

    def run(self, signals: Sequence[Signal2D]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x = self.get_x(signals)
        y = self.get_y(signals)
        z = np.zeros((len(signals), len(y), len(x)), dtype=signals[0].y.dtype)
        for i, sig in enumerate(signals):
            z[i, :, :] = ca_math.map_discrete_x_axis_2D(x, sig.x.astype(x.dtype), sig.z[:, :len(y)])
        return x, y, z
