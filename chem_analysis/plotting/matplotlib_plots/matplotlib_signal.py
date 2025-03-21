import numpy as np
import matplotlib.pyplot as plt

from chem_analysis.plotting.matplotlib_plots.matplotlib_utils import input_check
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.utils.math import get_slice
from chem_analysis.sec.sec_signal import SECSignal


def matplotlib_xy(
        fig: plt.Figure | None,
        plot_kwargs: dict,
        x: np.array,
        y: np.array,
        name: str,
        x_label: str,
        y_label: str,
        signal_: Signal
) -> plt.Figure:
    fig = input_check(fig)
    ax = fig.axes[0] if fig.axes else fig.add_subplot(111)

    kwargs = dict(linestyle="-", label=name)
    plot_kwargs = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    ax.plot(x, y, **plot_kwargs)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if isinstance(signal_, SECSignal):
        matplotlib_sec(x, y, signal_, ax)
    return fig


def matplotlib_sec(x: np.ndarray, y: np.ndarray, signal: SECSignal, ax: plt.axes):
    if signal.calibration is not None:
        bounds = signal.calibration.x_bounds
        if bounds[0] > bounds[1]:
            bounds = bounds[1], bounds[0]
        slice_ = get_slice(signal.x, bounds[0], bounds[1])
        max_ = np.max([2, np.max(y[slice_])])
        min_ = np.min([0, np.min(y[slice_])])
        span = (max_ - min_) * 0.05
        ax.ybound = [min_ - span, max_ + span]
        # ax.xbound = [0, 0.95]  # avoid overlap of legend and right y-axis


def matplotlib_discrete(
        fig: plt.Figure | None,
        plot_kwargs: dict,
        x: np.array,
        y: np.array,
        name: str,
        x_label: str,
        y_label: str,
        class_: type(Signal)
) -> plt.Figure:
    fig = input_check(fig)
    ax = fig.axes[0] if fig.axes else fig.add_subplot(111)

    kwargs = dict(x=x, y=y, name=name)
    plot_kwargs = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    ax.add_bar(**plot_kwargs)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    return fig
