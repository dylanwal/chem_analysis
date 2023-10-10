from typing import Callable

import numpy as np
import plotly.graph_objs as go  #TODO: change later

from chem_analysis.config import global_config
from chem_analysis.utils.plot_format import get_plot_color
from chem_analysis.sec.sec_signal import SECSignal
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.sec.sec_chromatogram import SECChromatogram
from chem_analysis.sec.sec_peak import SECPeak


class PlotSECConfig:
    def __init__(self):
        self.x_axis_mw_i = False  # False = x-axis is Time; True = x-axis is mw_i
        self.default_formatting: bool = True
        self.title = None
        self.y_label = None
        self.x_label = None

        # signal
        # self.signal_color: Callable | str = 'rgb(10,36,204)'
        self.signal_line_width = 2
        self._signal_color_counter = 0

        # peak
        # self.peak_show_trace: bool = False
        self.peak_show_shade: bool = True
        self.peak_show_bounds: bool = True
        self.peak_show_max: bool = True

        # calibration
        self.calibration_show_trace: bool = True
        self.calibration_show_limits: bool = True
        self.calibration_color = "rgb(130,130,130)"
        self.calibration_line_width = 1
        # self.calibration_shaded = False

    @property
    def peak_show(self) -> bool:
        vars_ = [getattr(self, k) for k in vars(self) if k.startswith("peak_show")]
        return any(vars_)

    @peak_show.setter
    def peak_show(self, value: bool):
        for k in vars(self):
            if k.startswith("peak_show"):
                setattr(self, k, value)

    @property
    def calibration_show(self) -> bool:
        vars_ = [getattr(self, k) for k in vars(self) if k.startswith("calibration_show")]
        return any(vars_)

    @calibration_show.setter
    def calibration_show(self, value: bool):
        for k in vars(self):
            if k.startswith("calibration_show"):
                setattr(self, k, value)

    def apply_default_formatting(self):
        ...

    def get_signal_color(self) -> str:
        return 'rgb(10,36,204)'

    def set_color_count(self):
        ...  # TODO


def plot_sec_chromatogram(chromatogram: SECChromatogram, config: PlotSECConfig = PlotSECConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_sec_chromatogram(chromatogram, config)

    raise NotImplementedError()


def plot_sec_signal(signal: SECSignal, config: PlotSECConfig = PlotSECConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        fig = go.Figure()
        plotly_sec_signal(fig, signal, config)
        plotly_sec_layout(fig, signal, config)
        return fig

    raise NotImplementedError()


def plotly_sec_chromatogram(chromatogram: SECChromatogram, config: PlotSECConfig):
    fig = go.Figure()
    if len(chromatogram.signals) == 0:
        return fig
    config.set_color_count(chromatogram.num_signals)

    for signal in chromatogram.signals:
        plotly_sec_signal(fig, signal, config)

    plotly_sec_layout(fig, signal, config)


def plotly_sec_layout(fig: go.Figure, signal: SECSignal, config: PlotSECConfig):
    if config.calibration_show and signal.calibration is not None:
        plotly_add_calibration(fig, signal, config)

    if config.default_formatting:
        config.apply_defualt_formatting(fig)

    # label axis
    if config.x_label is not None:
        fig.layout.xaxis.title = config.x_label
    else:
        if config.x_axis_mw_i:
            fig.layout.xaxis.title = "mw_i (g/mol)"
        else:
            fig.layout.xaxis.title = signal.x_label

    if config.y_label is not None:
        fig.layout.yaxis.title = config.y_label
    fig.layout.yaxis.title = signal.y_label


def plotly_sec_signal(fig: go.Figure, signal: SECSignal, config: PlotSECConfig) -> go.Figure:
    if config.peak_show:
        for peak in signal.peaks:
            plotly_add_peaks(fig, peak, config)

    # add main signal trace
    plot_kwargs = {
        "connectgaps": True,
        "legendgroup": "signal",
        "line": {
            "color": config.get_signal_color(),
            "width": config.signal_line_width
        }
    }
    fig.add_trace(
        go.Scatter(
            x=signal.mw_i if config.x_axis_mw_i else signal.x,
            y=signal.y,
            mode="lines",
            name=signal.name,
            **plot_kwargs
        )
    )

    return fig


def plotly_add_calibration(fig, signal: SECSignal, config: PlotSECConfig, legend_group: str = "calibration"):
    name = signal.calibration.name if signal.calibration.name is not None else "calibration"
    kwargs_ = {
        "name": name,
        "legendgroup": legend_group,
        "yaxis": "y2",
        "mode": "lines",
    }

    x = signal.mw_i if config.x_axis_mw_i else signal.x
    y = signal.calibration.get_y(x)

    # main calibration
    if config.calibration_show_trace:
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            line={"width": config.calibration_line_width, "color": config.calibration_color},
            showlegend=True,
            **kwargs_
        ))

    # add vertical limit bars
    if config.calibration_show_limits and not config.x_axis_mw_i:
        kwargs_["dash"] = 'dash'

        if config.x_axis_mw_i:
            x = (signal.calibration.lower_bound_y, signal.calibration.lower_bound_y)
        else:
            x = (signal.calibration.lower_bound_location, signal.calibration.lower_bound_location)
        fig.add_trace(go.Scatter(
            x=x,
            y=(0, np.max(y)),
            showlegend=False,
            line={"width": config.calibration_line_width, "color": config.calibration_color},
            **kwargs_
        ))

        if config.x_axis_mw_i:
            x = (signal.calibration.upper_bound_y, signal.calibration.upper_bound_y)
        else:
            x = (signal.calibration.upper_bound_location, signal.calibration.upper_bound_location)
        fig.add_trace(go.Scatter(
            x=x,
            y=(0, np.max(y)),
            showlegend=False,
            line={"width": config.calibration_line_width, "color": config.calibration_color},
            **kwargs_
        ))

    fig.update_layout(
        yaxis2=dict(
            title="molecular weight",
            titlefont=dict(
                color=config.calibration_color
            ),
            tickfont=dict(
                color=config.calibration_color
            ),
            anchor="x",
            overlaying="y",
            side="right",
            type="log",
            range=[2, 6]
        )
    )


def plotly_add_peaks(fig: go.Figure, peak: SECPeak, config: PlotSECConfig):
    if config.peak_show_shade:
        plotly_add_peak_shade(fig, peak, config)

    # if config.peak_show_trace:
    #     plotly_add_peak_trace()

    if config.peak_show_bounds:
        plotly_add_peak_bounds(fig, peak, config)

    if config.peak_show_max:
        plotly_add_peak_max(fig, peak, config)


def plotly_add_peak_shade(fig: go.Figure, peak, config):
    """ Plots the shaded area for the peak. """
    kwargs_ = {
        "width": 0
    }

    if group:
        kkwargs["legendgroup"] = group
    if y_label:
        kkwargs["yaxis"] = y_label

    fig.add_trace(go.Scatter(
        x=x,
        y=y,
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

def plot_mw(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True, y_label: str = None,
            title: str = None, spread: float = None, **kwargs) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    colors = get_plot_color(2)

    # add main trace
    plot_kwargs_xi = {
        "x": self.mw_i,
        "y": self.xi,
        "mode": 'lines',
        "connectgaps": True,
        "name": "xi",
        "line": dict(color=colors[0]),
        "yaxis": "y1"
    }

    fig.add_trace(go.Scatter(**plot_kwargs_xi))
    plot_kwargs_wi = {
        "x": self.mw_i,
        "y": self.wi,
        "mode": 'lines',
        "connectgaps": True,
        "name": "wi",
        "line": dict(color=colors[1]),
        "yaxis": "y2"
    }

    fig.add_trace(go.Scatter(**plot_kwargs_wi))

    # adding multiple y-axis
    y_axis_labels = ["y1", "y2"]
    if spread is None:
        spread = 0.05 * len(set(y_axis_labels))
    axis_format = get_multi_y_axis(colors, fig, spread)
    fig.update_layout(**axis_format)

    if auto_open:
        global fig_count
        fig.write_html(f'temp{fig_count}.html', auto_open=True)
        fig_count += 1

    return fig