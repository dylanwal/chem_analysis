import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.utils.math import get_slice
from chem_analysis.sec.sec_signal import SECSignal


def plotly_signal(
        signal: Signal,
        plot_kwargs: dict,
        fig: go.Figure | None,
        raw: bool = True,
        normalize: int = 0,
) -> go.Figure:
    fig = input_check(fig)

    if raw:
        name = signal.name + "_raw"
        x = signal.x_raw
        y = signal.y_raw
    else:
        x = signal.x
        y = signal.y
        name = signal.name

    if normalize == 1:
        y = signal.y_normalized_by_max()
    if normalize == 2:
        y = signal.y_normalized_by_area()

    if hasattr(signal, "_discrete"):
        plotly_signal_discrete_core(x, y, fig, name, plot_kwargs)
    else:
        plot_signal_core(x, y, fig, name, plot_kwargs)

    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    if isinstance(signal, SECSignal):
        plotly_signal_sec(signal, fig)
    return fig


def plot_signal_core(x: np.ndarray, y: np.ndarray, fig: go.Figure, name: str, plot_kwargs: dict):
    kwargs = dict(x=x, y=y, mode="lines", name=name)
    plot_kwargs = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    fig.add_scatter(**plot_kwargs)


def plotly_signal_sec(signal: SECSignal, fig: go.Figure):
    if signal.calibration is not None:
        bounds = signal.calibration.x_bounds
        if bounds[0] > bounds[1]:
            bounds = bounds[1], bounds[0]
        slice_ = get_slice(signal.x, bounds[0], bounds[1])
        max_ = np.max([2, np.max(signal.y[slice_])])
        min_ = np.min([0, np.min(signal.y[slice_])])
        span = (max_ - min_) * 0.05
        fig.layout.yaxis.range = [min_ - span, max_ + span]
        fig.layout.xaxis.domain = [0, 0.95]  # avoid overlap of legend and right y-axis


def plotly_signal_discrete_core(x: np.ndarray, y: np.ndarray, fig: go.Figure, name: str, plot_kwargs: dict):
    kwargs = dict(x=x, y=y, name=name)
    plot_kwargs = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    fig.add_bar(**plot_kwargs)
