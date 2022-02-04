from typing import Union, Callable
from enum import Enum

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from src.chem_analysis.analysis.base_obj.peak import Peak
from src.chem_analysis.analysis.algorithms import despike_methods, baseline_methods, smoothing_methods
import src.chem_analysis.analysis.algorithms as algorithms
from src.chem_analysis.analysis.logger import logger_analysis
from src.chem_analysis.analysis.utils import ObjList
from src.chem_analysis.analysis.utils.plot_format import get_plot_color, get_similar_color, add_plot_format

array_like = (np.ndarray, list, tuple)


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


class Signal:
    __count = 0

    def __init__(self,
                 name: str = None,
                 ser: pd.Series = None,
                 x: Union[np.ndarray, pd.Series] = None,
                 y: np.ndarray = None,
                 x_label: str = None,
                 y_label: str = None,
                 ):
        if name is None:
            name = f"trace_{Signal.__count}"
            Signal.__count += 1

        self.name = name
        self.x_label = x_label
        self.y_label = y_label

        if ser is not None and isinstance(ser, pd.Series) and x is None and y is None:
            self.raw = ser
            self.x_label = self.raw.index.name if self.x_label is not None else "x_axis"
            self.y_label = self.raw.index.name if self.y_label is not None else "y_axis"
        elif ser is None and isinstance(x, array_like) and isinstance(y, array_like):
            self.x_label = x_label if x_label is not None else "x_axis"
            self.y_label = y_label if y_label is not None else "y_axis"
            self.raw = pd.Series(data=y, index=x, name=self.y_label)
            self.raw.index.names = [self.x_label]

        if not hasattr(self, "raw"):
            mes = "Provide either a pandas Series (ser=) or two numpy arrays x and y (x=,y=)."
            raise ValueError(mes + f" (df:{type(ser)}, x:{type(x)}, y:{type(y)})")

        self.pipeline = ObjList(ProcessStep)
        self.peaks = ObjList(Peak)
        self._result = None
        self._result_up_to_date = False

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self.raw)})"
        return text

    @property
    def result(self) -> pd.Series:
        if not self._result_up_to_date:
            self.calc()
        return self._result

    @property
    def num_peaks(self) -> int:
        return len(self.peaks)

    def calc(self):
        x = self.raw.index.to_numpy()
        y = self.raw.to_numpy()
        for func in self.pipeline:
            x, y, z = func(x, y)

        self._result = pd.Series(y, x, name=self.y_label)
        self._result.index.names = [self.x_label]
        self._result_up_to_date = True

    def despike(self, method="default", **kwargs):
        """ despike

        Parameters
        ----------
        method
        kwargs

        Returns
        -------

        """
        if callable(method):
            func = method
        else:
            func: callable = despike_methods[method]

        self.pipeline.add(ProcessStep(func, ProcessType.DESPIKE, kwargs))
        self._result_up_to_date = False

        logger_analysis.debug(f"Despiking ({method}) done on: '{self.name}'.")

    def baseline(self, method="polynomial", **kwargs):
        if callable(method):
            func = method
        else:
            func: callable = baseline_methods[method]

        self.pipeline.add(ProcessStep(func, ProcessType.BASELINE, kwargs))
        self._result_up_to_date = False

        logger_analysis.debug(f"Baseline correction ({method}) done on: '{self.name}'.")

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

        logger_analysis.debug(f"Smooth ({method}) done on: '{self.name}'.")

    def auto_peak_picking(self, **kwargs):
        kwargs_ = {"width": self.raw.index[-1] / 50}
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}
        peaks_index = algorithms.scipy_find_peaks(self.result.to_numpy(), **kwargs_)
        if len(peaks_index) != 0:
            for peak in peaks_index:
                lb, ub = algorithms.rolling_value(self.result.to_numpy(), peak_index=peak, sensitivity=0.5,
                                                  cut_off=0.05)
                self.peaks.add(Peak(self, lb, ub))
        else:
            logger_analysis.warning(f"No peaks found in signal '{self.name}'.")

        logger_analysis.debug(f"Auto peak picking done on: '{self.name}'. Peaks found: {self.num_peaks}")

    def stats(self, op_print: bool = True):
        text = ""
        for i, peak in enumerate(self.peaks):
            if i == 0:
                text += peak.stats(op_print=False)
            else:
                text += peak.stats(op_print=False, op_headers=False)

        if op_print:
            print(text)
        return text

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, **kwargs) -> go.Figure:
        """ Basic plotting """
        if fig is None:
            fig = go.Figure()

        if "color" in kwargs:
            color = kwargs.pop("color")
        else:
            color = 'rgb(0,0,0)'

        group = self.name

        # add peaks
        if op_peaks:
            if len(self.peaks) > 0:
                if color == 'rgb(0,0,0)':
                    peak_color = get_plot_color(self.num_peaks)
                else:
                    peak_color = get_similar_color(color, self.num_peaks)

                for peak, color_ in zip(self.peaks, peak_color):
                    peak._plot(fig, color=color_, group=group)

        # add main trace
        fig.add_trace(
            go.Scatter(x=self.result.index, y=self.result, mode="lines", connectgaps=True, name=self.result.name,
                       line=dict(color=color), legendgroup=group))

        if auto_format:
            add_plot_format(fig, self.result.index.name, self.result.name)

        if auto_open:
            fig.write_html(f'temp.html', auto_open=True)

        return fig


def local_run():
    from scipy.stats import norm
    n = 1000
    rv = norm(loc=n/2, scale=10)
    x = np.linspace(0, n, n)
    y = np.linspace(0, n, n) + 20 * np.random.random(n) + 5000 * rv.pdf(x)
    trace = Signal(name="test", x=x, y=y, x_label="x_test", y_label="y_test")
    trace.baseline(deg=1)
    trace.auto_peak_picking()
    trace.stats()
    trace.plot()
    print("done")


if __name__ == '__main__':
    local_run()
