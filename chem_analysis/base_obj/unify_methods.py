import abc
from typing import Sequence

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import chem_analysis.utils.math as ca_math
from chem_analysis.base_obj.signal_ import Signal


class UnifyMethod(abc.ABC):
    @abc.abstractmethod
    def run(self, signals: Sequence[Signal]) -> tuple[np.ndarray, np.ndarray]:
        """ return x (shape n) and z (shape m x n) """
        ...


class UnifyMethodStrict(UnifyMethod):
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

    def run(self, signals: Sequence[Signal]) -> tuple[np.ndarray, np.ndarray]:
        x = signals[0].x
        z = np.empty((len(signals), len(x)), dtype=signals[0].y.dtype)
        args = self.get_args()
        for i, sig in enumerate(signals):
            if not np.all(np.isclose(sig.x, x, **args)):
                raise ValueError(f"Signal {i} has a different x-axis than first signal.")
            z[i, :] = sig.y

        return x, z


class UnifyMethodExpandValueStrict(UnifyMethod):
    """
    Uses the smallest min and largest max
    fills with a value
    """
    def __init__(self,
                 value: int | float = 0,
                 min_: None | int | float = None,
                 max_: None | int | float = None,
                 rtol: float | None = None,
                 atol: float | None = None
                 ):
        """

        Parameters
        ----------
        value: int | float
            value to fill in expansion with
        min_:
            set to override the smallest min
        max_:
            set to override the largest max
        """
        self.min_ = min_
        self.max_ = max_
        self.value = value
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

    def get_x(self, signals: Sequence[Signal]) -> np.ndarray:
        if self.max_ is None:
            max_ = np.max([sig.x[-1] for sig in signals])
        else:
            max_ = self.max_
        if self.min_ is None:
            min_ = np.min([sig.x[0] for sig in signals])
        else:
            min_ = self.min_

        # try to get an existing x-axis
        for sig in signals:
            if sig.x[0] == min_ and sig.x[1] == max_:
                return sig.x

        # try to compose existing x-axis
        return ca_math.get_all_unique_values_within_tolerance([sig.x for sig in signals], **self.get_args())

    def run(self, signals: Sequence[Signal]) -> tuple[np.ndarray, np.ndarray]:
        x = self.get_x(signals)

        z = np.ones((len(signals), len(x)), dtype=signals[0].y.dtype)*self.value
        args = self.get_args()
        for i, sig in enumerate(signals):
            if np.all(np.isclose(sig.x, x, **args)):
                z[i, :] = sig.y
            else:
                # loop through to find where sig.x aligns with x
                offset = 0
                for ii, x_ in enumerate(x):
                    if np.isclose(sig.x[0], x_, **args):
                        offset = ii
                        break
                    else:
                        raise ValueError(f"Signal {i} has a different x-axis than first signal.")
                z[i, offset:offset+len(sig.y)] = sig.y

        return x, z


class UnifyMethodExpandInterpolate(UnifyMethod):
    """
    Uses the smallest min and largest max
    fills with a value
    """
    def __init__(self,
                 value: None | int | float = None,
                 min_: None | int | float = None,
                 max_: None | int | float = None
                 ):
        """

        Parameters
        ----------
        value: None | int | float
            None will extrapolate outside bounds
            value to fill in expansion with
        min_:
            set to override the smallest min
        max_:
            set to override the largest max
        """
        self.min_ = min_
        self.max_ = max_
        self.value = value

    def get_x(self, signals: Sequence[Signal]) -> np.ndarray:
        if self.max_ is None:
            max_ = np.max([sig.x[-1] for sig in signals])
        else:
            max_ = self.max_
        if self.min_ is None:
            min_ = np.min([sig.x[0] for sig in signals])
        else:
            min_ = self.min_

        # try to get an existing x-axis
        for sig in signals:
            if sig.x[0] == min_ and sig.x[1] == max_:
                return sig.x

        # create new x-axis
        max_data_points = np.max([len(sig.x) for sig in signals])
        return np.linspace(min_, max_, max_data_points)

    def run(self, signals: Sequence[Signal]) -> tuple[np.ndarray, np.ndarray]:
        x = self.get_x(signals)

        z = np.ones((len(signals), len(x)), dtype=signals[0].y.dtype)*self.value
        for i, sig in enumerate(signals):
            if not np.all(np.isclose(sig.x, x, rtol=0.0001)):
                z[i, :] = sig.y
            else:
                spline = InterpolatedUnivariateSpline(x=sig.x, y=sig.y)
                z[i, :] = spline(x)

        return x, z


class MethodShrink(UnifyMethod):
    """
        Uses the largest min and smallest max
        cuts values
    """
    def __init__(self, min_: None | int | float = None, max_: None | int | float = None):
        """

        Parameters
        ----------
        min_:
            set to override the largest min
        max_:
            set to override the smallest max
        """
        self.min_ = min_
        self.max_ = max_

    def run(self, data: Sequence[np.ndarray]) -> tuple[int | float, int | float]:
        """

        Parameters
        ----------
        data:
            data must be sorted

        Returns
        -------

        """
        if self.max_ is None:
            max_ = np.min([sig[-1] for sig in data])
        else:
            max_ = self.max_
        if self.min_ is None:
            min_ = np.max([sig[0] for sig in data])
        else:
            min_ = self.min_

        return min_, max_
