from typing import Union, Callable
from enum import Enum
from functools import wraps

import numpy as np
import pandas as pd

from src.chem_analysis.analysis.algorithms import despike_methods, baseline_methods, smoothing_methods, \
    peak_picking_methods, bound_detection_methods
from src.chem_analysis.analysis.utils import ObjList, Pipeline
from src.chem_analysis.analysis.utils.plot_format import plot_series


array_like = (np.ndarray, list, tuple)


class Peak:
    def __init__(self, id_: int, low_bound: float, high_bound: float):
        self.id_ = id_

        self.low_bound = low_bound
        self.high_bound = high_bound

        self.max = None
        self.area = None
        self.std = None
        self.skew = None
        self.k = None

    def __repr__(self):
        return f"peak: {self.id_} at {self.max}"


class ProcessType(Enum):
    RAW = 0
    DESPIKE = 1
    BASELINE = 2
    SMOOTH = 3


class ProcessStep:
    def __init__(self, func: Callable, type_: ProcessType, kwargs: dict = None):
        self.func = func
        self.current_type = type_
        self.kwargs = kwargs

    def __repr__(self):
        return f"{self.func.__name__}; type:{self.current_type}"

    def __call__(self, *args, **kwargs):
        if kwargs is not None:
            kwargs = {**kwargs, **self.kwargs}

        return self.func(*args, **kwargs)


class Trace:
    despike_methods = despike_methods
    baseline_methods = baseline_methods
    smoothing_methods = smoothing_methods

    def __init__(self,
                 ser: pd.Series = None,
                 x: Union[np.ndarray, pd.Series] = None,
                 y: np.ndarray = None,
                 x_label: str = None,
                 y_label: str = None,
                 ):

        self.x_label = x_label
        self.y_label = y_label

        if ser is not None and isinstance(ser, pd.Series) and x is None and y is None:
            self.ser = ser
            self.x_label = self.ser.index.name if self.x_label is None else "x_axis"
            self.y_label = self.ser.index.name if self.y_label is None else "y_axis"
        elif ser is None and isinstance(x, array_like) and isinstance(y, array_like):
            self.x_label = x_label if x_label is not None else "x_axis"
            self.y_label = y_label if y_label is not None else "y_axis"
            self.ser = pd.Series(data=y, index=x, name=self.y_label)
            self.ser.index.names = [self.x_label]

        if not hasattr(self, "ser"):
            mes = "Provide either a pandas Series (df=) or two numpy arrays x and y (x=,y=)."
            raise ValueError(mes + f" (df:{type(ser)}, x:{type(x)}, y:{type(y)})")

        self.pipeline = Pipeline()
        self.peaks = ObjList(Peak)
        self._result = None
        self._result_up_to_date = False

    def __repr__(self):
        text = ""
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self.ser)})"
        return text

    @property
    def result(self):
        if not self._result_up_to_date:
            self.calc()
        return self._result

    def calc(self):
        x = self.ser.index.to_numpy()
        y = self.ser.to_numpy()
        for func in self.pipeline:
            x, y, z = func(x, y)

        self._result = pd.Series(y, x, name=self.y_label)
        self._result.index.names = [self.x_label]
        self._result_up_to_date = True

    def despike(self, method="default", **kwargs):
        if callable(method):
            func = method
        else:
            func: callable = despike_methods[method]

        self.pipeline.add(ProcessStep(func, ProcessType.DESPIKE, kwargs))
        self._result_up_to_date = False

    def baseline(self, method="polynomial", **kwargs):
        if callable(method):
            func = method
        else:
            func: callable = baseline_methods[method]

        self.pipeline.add(ProcessStep(func, ProcessType.BASELINE, kwargs))
        self._result_up_to_date = False

    def smooth(self, method="default", **kwargs):
        if callable(method):
            func = method
        else:
            try:
                func: callable = smoothing_methods[method]
            except KeyError:
                raise ValueError(f"Not a valid 'peak_picking' method. (given: {method})")

        self.pipeline.add(ProcessStep(func, ProcessType.SMOOTH, kwargs))
        self._result_up_to_date = False

    def peak_picking(self, method_peak="scipy", method_bound="rolling_value", **kwargs):
        if callable(method_peak):
            func_peak = method_peak
        else:
            try:
                func_peak: callable = peak_picking_methods[method_peak]
            except KeyError:
                raise ValueError(f"Not a valid 'peak_picking' method. (given: {method_peak})")

        if callable(method_bound):
            func_bound = method_bound
        else:
            try:
                func_bound: callable = bound_detection_methods[method_bound]
            except KeyError:
                raise ValueError(f"Not a valid 'bound_detection' method. (given: {method_bound})")

        peaks = func_peak(self.result, **kwargs)
        bounds = []
        for peak in peaks:
            lb, ub = func_bound(self.result, peak_index=peak,)


    @wraps(plot_series)
    def plot(self, **kwargs):
        """ Basic plotting """
        return plot_series(data=self.result, **kwargs)


def local_run():
    trace = Trace(
        x=np.linspace(0, 100, 100),
        y=np.linspace(0, 100, 100) + 30 * np.random.random(100) + 50 * np.random.normal(50, 10, 100),
        x_label="x_test",
        y_label="y_test"
    )
    trace.baseline(deg=1)
    trace.plot()
    print("done")


if __name__ == '__main__':
    local_run()
