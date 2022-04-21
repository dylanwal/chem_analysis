from typing import Protocol, Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go

# from chem_analysis.analysis.utils.plot_format import add_plot_format
from chem_analysis.analysis.utils.sig_fig import sig_figs
from chem_analysis.analysis.logger import logger_analysis
from chem_analysis.analysis.utils import fig_count


class PeakSupports(Protocol):
    name: str
    result: pd.Series


def _length_limit(label: str, limit: int) -> str:
    """ For printing out stats. """
    if len(label) > limit:
        return label[0:limit - 5] + "..."

    return label


class Peak:
    """ Peak

    a peak of a signal

    Attributes
    ----------
    id_: int
        id of peak
    x: nd.ndarray
        x axis data
    y: nd.ndarray
        y axis data
    slice: slice
        slice of full signal for peak
    lb_index: int
        index of lower bound
    hb_index: int
        index of higher bound
    lb_loc: float
        x location of lower bound
    hb_loc: float
        x location of higher bound
    lb_value: float
        y value of lower bound
    hb_value: float
        y value of higher bound
    max_index: int
        index of peak maximum
    max: float
        y value of peak maximum
    max_loc: float
        x location of peak maximum
    area: float
        area under the peak
    mean: float
        average value
    std: float
        standard deviation
    skew: float
        skew
        symmetric: -0.5 to 0.5; moderate skew: -1 to -0.5 or 0.5 to 1; high skew: <-1 or >1;
        positive tailing to higher numbers; negative tailing to smaller numbers
    kurtosis: float
        kurtosis (Fisher) (Warning: highly sensitive to peak bounds)
        negative: flatter peak; positive: sharp peak
    fwhm: float
        full width at half maximum
    asym: float
        asymmetry factor; distance from the center line of the peak to the back slope divided by the distance from the
        center line of the peak to the front slope;
        >1 tailing to larger values; <1 tailing to smaller numbers

    """
    def __init__(self, parent: PeakSupports, lb_index: int, hb_index: int, id_: int = None):
        self.id_ = id_
        self._parent = parent

        self.x = None
        self.y = None
        self.y_norm = None

        self._lb_index = None
        self._hb_index = None
        self.lb_loc = None
        self.hb_loc = None
        self.lb_value = None
        self.hb_value = None
        self.max_index = None
        self.max = None
        self.max_loc = None
        self.area = None

        self.mean = None
        self.std = None
        self.skew = None
        self.kurtosis = None
        self.fwhm = None
        self.asym = None

        self.lb_index = lb_index
        self.hb_index = hb_index

    def __repr__(self):
        return f"peak: {self.id_} at {self.max_loc}"

    @property
    def lb_index(self):
        return self._lb_index

    @lb_index.setter
    def lb_index(self, lb_index):
        self._lb_index = lb_index
        if self._hb_index is not None:
            self.calc()

    @property
    def hb_index(self):
        return self._hb_index

    @hb_index.setter
    def hb_index(self, hb_index):
        self._hb_index = hb_index
        if self._lb_index is not None:
            self.calc()

    @property
    def slice(self):
        return slice(self.lb_index, self.hb_index)

    def calc(self):
        """ Calculates the stats for the peak. """
        self.x = self._parent.result.index[self.slice].to_numpy()
        self.y = self._parent.result.iloc[self.slice].to_numpy()
        self.y_norm = self.y/np.trapz(x=self.x, y=self.y)
        self.lb_value = self.y[0]
        self.hb_value = self.y[-1]
        self.lb_loc = self.x[0]
        self.hb_loc = self.x[-1]

        self.max_index = np.argmax(self.y) + self.lb_index
        self.max = self.y[self.max_index-self.lb_index]
        self.max_loc = self.x[self.max_index-self.lb_index]

        self.mean = np.trapz(x=self.x, y=self.x*self.y_norm)
        self.std = np.sqrt(np.trapz(x=self.x, y=self.y_norm*(self.x-self.mean)**2))
        self.skew = np.trapz(x=self.x, y=self.y_norm * (self.x - self.mean) ** 3) / self.std ** 3
        self.kurtosis = (np.trapz(x=self.x, y=self.y_norm * (self.x - self.mean) ** 4) / self.std ** 4) - 3
        self.fwhm = self.get_fw(x=self.x, y=self.y, height=0.5)
        self.asym = self.get_asym(x=self.x, y=self.y, height=0.1)

        self.area = np.trapz(x=self.x, y=self.y)

    def get_fw(self, x: np.ndarray, y: np.ndarray, height: float = 0.5) -> Union[float, int]:
        """ Calculates full width at a height. """
        lower, high = self.get_width_at(x, y, height)
        return abs(high-lower)

    def get_asym(self, x: np.ndarray, y: np.ndarray, height: float = 0.1) -> Union[float, int]:
        """ Calculates asymmetry factor at height. """
        lower, high = self.get_width_at(x, y, height)
        middle = x[np.argmax(y)]

        return (high-middle)/(middle-lower)

    @staticmethod
    def get_width_at(x: np.ndarray, y: np.ndarray, height: float = 0.5) -> tuple[float, float]:
        """ Determine full-width-x_max of a peaked set of points, x and y. """
        height_half_max = np.max(y) * height
        index_max = np.argmax(y)
        if index_max == 0 or index_max == len(x):  # peak max is at end.
            logger_analysis.info("Finding fwhm is not possible with a peak max at an bound.")
            return 0, 0

        x_low = np.interp(height_half_max, y[:index_max], x[:index_max])
        x_high = np.interp(height_half_max, np.flip(y[index_max:]), np.flip(x[index_max:]))

        if x_low == x[0]:
            logger_analysis.info("fwhm or asym is having to linear interpolate on the lower end.")
            slice_ = max(3, int(index_max/10))
            fit = np.polyfit(y[:slice_], x[:slice_], deg=1)
            p = np.poly1d(fit)
            x_low = p(height_half_max)

        if x_high == x[-1]:
            logger_analysis.info("fwhm or asym is having to linear interpolate on the lower end.")
            slice_ = max(3, int(index_max/10))
            fit = np.polyfit(y[-slice_:], x[-slice_:], deg=1)
            p = np.poly1d(fit)
            x_high = p(height_half_max)

        return x_low, x_high

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             y_label: str = None, title: str = None, **kwargs) -> go.Figure:
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

        # add main trace
        plot_kwargs = {
            "x": self.x,
            "y": self.y,
            "mode": 'lines',
            "connectgaps": True,
            "name": self._parent.result.name,
            "line": dict(color=color)
        }
        if y_label is not None:
            plot_kwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(**plot_kwargs))

        if auto_format:
            if title is not None:
                fig.update_layout(title=title)
            add_plot_format(fig, self._parent.result.index.name, str(self._parent.result.name))

        if auto_open:
            global fig_count
            fig.write_html(f'temp{fig_count}.html', auto_open=True)
            fig_count += 1

        return fig

    def plot_add_on(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plot

        Plots several things related to a peak.

        Parameters
        ----------
        fig: go.Figure
            plotly figure
        group: str
            the legend table group that the plot_add_on will be associated with (need to make things appear/disappear together)
        y_label: str
            Label for the y-axis the data is associated with

        """
        self._plot_max(fig, group, y_label, **kwargs)
        # self._plot_bounds(fig, group, **kwargs)
        self._plot_shade(fig, group, y_label, **kwargs)

    def _plot_shade(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plots the shaded area for the peak. """
        kwargs_ = {
            "width": 0
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(
            x=self.x,
            y=self.y,
            mode="lines",
            fill='tozeroy',
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))

    def _plot_max(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plots peak name at max. """
        kwargs_ = {
            "size": 1
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(
            x=[self.max_loc],
            y=[self.max],
            mode="text",
            marker=kwargs_,
            text=[f"{self._parent.name}: {self.id_}"],
            textposition="top center",
            showlegend=False,
            **kkwargs
        ))

    def _plot_bounds(self, fig: go.Figure, group: str = None, height: float = 0.08, y_label: str = None, **kwargs):
        """ Adds bounds at the bottom of the plot_add_on for peak area. """
        kwargs_ = {
            "width": 5
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        bound_height = max(self._parent.result) * height
        # bounds
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.lb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs

        ))
        fig.add_trace(go.Scatter(
            x=[self.hb_loc, self.hb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.hb_loc],
            y=[0, 0],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))

    def stats(self, op_print: bool = True, op_headers: bool = True, window: int = 150, headers: dict = None):
        """ Prints stats out for peak. """
        text = ""
        if headers is None:
            headers = {  # attribute: print
                "id_": "id", "lb_loc": "low bound", "max_loc": "max", "hb_loc": "high bound", "area": "area"
            }

        # format
        width = int(window / len(headers))
        row_format = ("{:<" + str(width) + "}") * len(headers)

        # headers
        if op_headers:
            headers_ = [_length_limit(head, width) for head in headers.values()]
            text = row_format.format(*headers_) + "\n"
            text = text + "-" * window + "\n"

        # values
        entries = [_length_limit(str(sig_figs(getattr(self, k), 3)), width) for k in headers]
        entries[0] = f"{self._parent.name}: {entries[0]}"  # peak name: peak id
        text = text + row_format.format(*entries) + "\n"

        if op_print:
            print(text)

        return text
