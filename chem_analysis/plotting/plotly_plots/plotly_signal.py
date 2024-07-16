import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_discrete import SignalDiscrete
from chem_analysis.utils.math import get_slice
from chem_analysis.sec.sec_signal import SECSignal


def plotly_signal(
        signal: Signal | SignalDiscrete,
        fig: go.Figure | None,
        config: PlotlyConfig | None = None,
        raw: bool = True
) -> go.Figure:
    fig, config = PlotlyConfig.input_check(fig, config)

    if raw:
        name = signal.name + "_raw"
        x = signal.x_raw
        y = signal.y_raw
    else:
        x = signal.x
        y = signal.y
        name = signal.name

    if config.normalize is config.NORMALIZATION_OPTIONS.PEAK_HEIGHT:
        y = signal.y_normalized_by_max()

    if isinstance(signal, SignalDiscrete):
        plotly_signal_discrete_core(x, y, fig, name)
    else:
        plot_signal_core(x, y, fig, name)

    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    if isinstance(signal, SECSignal):
        plotly_signal_sec(signal, fig)
    return fig


def plot_signal_core(x: np.ndarray, y: np.ndarray, fig: go.Figure, name: str):
    fig.add_scatter(
        x=x,
        y=y,
        mode="lines",
        name=name,
    )


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


def plotly_signal_discrete_core(x: np.ndarray, y: np.ndarray, fig: go.Figure, name: str):
    fig.add_bar(
        x=x,
        y=y,
        name=name,
    )
