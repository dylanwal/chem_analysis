from typing import Union, Callable
from enum import Enum

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.base_obj.peak import Peak
from chem_analysis.analysis.algorithms import despike_methods, baseline_methods, smoothing_methods
import chem_analysis.analysis.algorithms as algorithms
from chem_analysis.analysis.logger import logger_analysis
from chem_analysis.analysis.utils import ObjList, fig_count
from chem_analysis.analysis.utils.plot_format import get_plot_color, get_similar_color, add_plot_format

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
    """ signal

    A signal is any x-y data.

    Attributes
    ----------
    name: str
        Any name the user wants to add.
    raw: pd.Series
        raw data
    x_label: str
        x-axis label
    y_label: str
        y-axis label
    pipeline:
        data processing pipeline
    peaks: Peak
        peaks found in signal
    num_peaks: int
        number of peaks
    result: pd.Series
        data process through pipeline
    result_norm: pd.Series
        data process through pipeline normalized (peak max = 1)

    """
    __count = 0
    _peak = Peak

    def __init__(self,
                 name: str = None,
                 ser: pd.Series = None,
                 x: Union[np.ndarray, pd.Series] = None,
                 y: np.ndarray = None,
                 x_label: str = None,
                 y_label: str = None,
                 _parent=None
                 ):
        """

        Parameters
        ----------
        name: str
            user defined name
        ser: pd.Series
            x-y data in the form of a pandas Series
        x: np.ndarray
            x data
        y: np.ndarray
            y data
        x_label: str
            x-axis label
        y_label: str
            y-axis label

        Notes
        -----
        * Either 'ser' or 'x' and 'y' are required but not both.

        """
        if name is None:
            name = f"trace_{Signal.__count}"
            Signal.__count += 1

        self.name = name
        self.x_label = x_label
        self.y_label = y_label

        if ser is not None and isinstance(ser, pd.Series) and x is None and y is None:
            self.raw = ser
            if self.x_label is None:
                self.x_label = self.raw.index.name
            if self.y_label is None:
                self.y_label = self.raw.name
        elif ser is None and isinstance(x, array_like) and isinstance(y, array_like):
            self.raw = pd.Series(data=y, index=x, name=self.y_label)
            self.raw.index.names = [self.x_label]

        if not hasattr(self, "raw"):
            mes = "Provide either a pandas Series (ser=) or two numpy arrays x and y (x=,y=)."
            raise ValueError(mes + f" (df:{type(ser)}, x:{type(x)}, y:{type(y)})")

        if self.x_label is None:
            self.x_label = "x_axis"
        if self.y_label is None:
            self.y_label = "y_axis"

        self.pipeline = ObjList(ProcessStep)
        self.peaks = ObjList(self._peak)
        self._result = None
        self._result_norm = None
        self._result_up_to_date = False
        self._parent = _parent

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
    def result_norm(self) -> pd.Series:
        if not self._result_up_to_date:
            self.calc()
        return self._result_norm

    @property
    def num_peaks(self) -> int:
        return len(self.peaks)

    def calc(self):
        x = self.raw.index.to_numpy()
        y = self.raw.to_numpy()
        for func in self.pipeline:
            x, y, z = func(x, y)

        self._result = pd.Series(y, x, name=self.y_label)
        self._result_norm = self._result/np.max(self._result)
        self._result.index.names = [self.x_label]
        self._result_up_to_date = True

    def despike(self, method="default", **kwargs):
        if callable(method):
            func = method
        else:
            try:
                func: callable = despike_methods[method]
            except KeyError:
                raise ValueError(f"Not a valid 'despiking' method. (given: {method})")

        self.pipeline.add(ProcessStep(func, ProcessType.DESPIKE, kwargs))
        self._result_up_to_date = False

        logger_analysis.debug(f"Despiking ({method}) done on: '{self.name}'.")

    def baseline(self, method="polynomial", **kwargs):
        if callable(method):
            func = method
        else:
            try:
                func: callable = baseline_methods[method]
            except KeyError:
                raise ValueError(f"Not a valid 'baseline' method. (given: {method})")

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
                raise ValueError(f"Not a valid 'smooth' method. (given: {method})")

        self.pipeline.add(ProcessStep(func, ProcessType.SMOOTH, kwargs))
        self._result_up_to_date = False

        logger_analysis.debug(f"Smooth ({method}) done on: '{self.name}'.")

    def auto_peak_picking(self, limit_range: list[float] = None, **kwargs):
        self.peaks.clear()

        # find peaks
        kwargs_ = {"width": self.raw.index[-1] / 100, "height": 0.03, "prominence": 0.03}
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}
        if limit_range:
            lb_index = np.argmin(np.abs(self.result.index.to_numpy() - limit_range[0]))
            ub_index = np.argmin(np.abs(self.result.index.to_numpy() - limit_range[1]))
            y = self.result_norm.iloc[lb_index:ub_index].to_numpy()
            y = y/np.max(y)
            peaks_index = algorithms.scipy_find_peaks(y, **kwargs_) + lb_index
        else:
            peaks_index = algorithms.scipy_find_peaks(self.result_norm.to_numpy(), **kwargs_)

        # get bounds from peak maximums
        if len(peaks_index) != 0:
            for peak in peaks_index:
                lb, ub = algorithms.rolling_value(self.result_norm.to_numpy(), peak_index=peak, sensitivity=0.1,
                                                  cut_off=0.05)
                self.peaks.add(self._peak(self, lb, ub))
        else:
            logger_analysis.warning(f"No peaks found in signal '{self.name}'.")

        logger_analysis.debug(f"Auto peak picking done on: '{self.name}'. Peaks found: {self.num_peaks}")

    def auto_peak_baseline(self, iterations: int = 3, limit_range: list[float] = None, **kwargs):
        """
        Does automatic baseline correction and peak detection.
        It uses peak detection to create a mask for baseline correction in an iterative fashion.

        """
        self.baseline(**kwargs)
        self.auto_peak_picking(limit_range=limit_range)
        for i in range(iterations):
            self.pipeline.remove(-1)
            mask = np.ones_like(self.raw.to_numpy())
            for peak in self.peaks:
                mask[peak.slice] = False
            self.baseline(mask=mask, **kwargs)
            self.auto_peak_picking()
            logger_analysis.debug(f"'auto_peak_baseline' iteration {i} done.")

    def stats(self, op_print: bool = True, op_headers: bool = True) -> str:
        """ Print out signal/peak stats. """
        text = ""
        for i, peak in enumerate(self.peaks):
            if i == 0:
                if op_headers:
                    text += peak.stats(op_print=False)
                    continue

            text += peak.stats(op_print=False, op_headers=False)

        if op_print:
            print(text)
        return text

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, y_label: str = None, title: str = None, **kwargs) -> go.Figure:
        """ Plot

        General plotting

        Parameters
        ----------
        fig: go.Figure
            plotly figure; will automatically create if not provided
        auto_open: bool
            create "temp.html" and auto_open in browser
        auto_format: bool
            apply built-in formatting
        op_peaks: bool
            add peak plotting stuff
        y_label: str
            y_axis label (used for multiple y-axis)
        title: str
            title

        Returns
        -------
        fig: go.Figure
            plotly figure

        """
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
                    peak.plot_add_on(fig, color=color_, group=group, y_label=y_label)

        # add main trace
        plot_kwargs = {
            "x": self.result.index,
            "y": self.result,
            "mode": 'lines',
            "connectgaps": True,
            "name": self.result.name,
            "legendgroup": group,
            "line": dict(color=color)
        }
        if y_label is not None:
            plot_kwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(**plot_kwargs))

        if auto_format:
            if title is not None:
                fig.update_layout(title=title)
            add_plot_format(fig, self.result.index.name, str(self.result.name))

        if auto_open:
            global fig_count
            fig.write_html(f'temp{fig_count}.html', auto_open=True)
            fig_count += 1

        return fig


def local_run():
    from scipy.stats import norm
    n = 1000
    rv = norm(loc=n/2, scale=10)
    x = np.linspace(0, n, n)
    y = np.linspace(0, n, n) + 20 * np.random.random(n) + 5000 * rv.pdf(x)
    signal = Signal(name="test", x=x, y=y, x_label="x_test", y_label="y_test")
    signal.baseline(deg=1)
    signal.auto_peak_picking()
    signal.stats()
    signal.plot()
    print("done")


if __name__ == '__main__':
    local_run()
