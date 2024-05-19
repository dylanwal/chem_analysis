
import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.utils.math import get_slice
from chem_analysis.sec.sec_signal import SECSignal


def plotly_signal(signal: Signal, fig: go.Figure | None, config: PlotlyConfig | None = None, raw: bool = True) \
        -> go.Figure:
    fig, config = PlotlyConfig.input_check(fig, config)

    # can add if statements for special plotting per signal
    # if isinstance(signal, SECSignal):
    #     return

    if raw:
        name = signal.name + "_raw"
        x = signal.x_raw
        y = signal.data_raw
    else:
        x = signal.x
        y = signal.y
        name = signal.name

    if config.normalize is config.NORMALIZATION_OPTIONS.PEAK_HEIGHT:
        y = signal.y_normalized_by_max()

    plot_signal_core(x, y, fig, name)
    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    if isinstance(signal, SECSignal):
        plotly_signal_sec(signal, fig)
    return fig


def plot_signal_core(x: np.ndarray, y: np.ndarray, fig:go.Figure, name: str):
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=name,
        )
    )


def plotly_signal_sec(signal: SECSignal, fig: go.Figure):
    if signal.calibration is not None:
        slice_ = get_slice(signal.x, *signal.calibration.x_bounds)
        max_ = np.max([2, np.max(signal.y[slice_])])
        min_ = np.min([0, np.min(signal.y[slice_])])
        span = (max_ - min_) * 0.05
        fig.layout.yaxis.range = [min_-span, max_+span]
