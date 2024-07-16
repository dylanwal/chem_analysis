import abc
from typing import Sequence

import numpy as np

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
