import copy
import logging
from typing import Sequence

import numpy as np

from chem_analysis.config import global_config
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.processing.processor import Processor
from chem_analysis.processing.processing_method import Baseline

logger = logging.getLogger(__name__)


def signal(
        signal_: Signal,
        *,
        fig=None,
        normalize: int = 0,
        plot_kwargs: dict | None = None,
):
    """

    Parameters
    ----------
    signal_:
        signal to be plotted
    fig:
        figure for signal to be added to
    normalize:
        0: no normalization
        1: normalize by height
        2: normalize by area
    plot_kwargs:
        arg passed directly to plotting function

    Returns
    -------

    """
    plotting_lib = global_config.get_plotting_lib() if fig is None else global_config.get_plotting_lib_from_fig(fig)

    x = signal_.x
    if normalize == 1:
        y = signal_.y_normalized_by_max()
    elif normalize == 2:
        y = signal_.y_normalized_by_area()
    else:
        y = signal_.y

    name = signal_.name
    x_label = signal_.x_label
    y_label = signal_.y_label

    discrete_flag = False
    if hasattr(signal_, "_discrete"):
        discrete_flag = True

    plot_kwargs = copy.copy(plot_kwargs) or {}
    if plotting_lib == global_config.PLOTTING_LIBRARIES.plotly:
        from chem_analysis.plotting.plotly_plots.plotly_signal import plotly_xy, plotly_discrete
        if discrete_flag:
            return plotly_discrete(fig, plot_kwargs, x, y, name, x_label, y_label, signal_)
        else:
            return plotly_xy(fig, plot_kwargs, x, y, name, x_label, y_label, signal_)

    if plotting_lib == global_config.PLOTTING_LIBRARIES.matplotlib:
        from chem_analysis.plotting.matplotlib_plots.matplotlib_signal import matplotlib_xy, matplotlib_discrete
        if discrete_flag:
            return matplotlib_discrete(fig, plot_kwargs, x, y, name, x_label, y_label, signal_)
        else:
            return matplotlib_xy(fig, plot_kwargs, x, y, name, x_label, y_label, signal_)

    if plotting_lib == global_config.PLOTTING_LIBRARIES.pygraphqt:
        raise NotImplementedError("not done yet")

    raise RuntimeError('Should never happen')


def add_peaks(
        signal_: Signal,
        bounds: np.ndarray | None = None,
        peaks: np.ndarray | None = None,
        *,
        fig=None,
        labels: Sequence[str] = None,
        label_mode: int | Sequence[int] = 0,
        mode: int | Sequence[int] = 0,
        normalize: int = 0,
        plot_kwargs: dict | None = None,
):
    """

    Parameters
    ----------
    signal_
    bounds:
        a [n,2] array with index of [[left, right], [left, right], ...]
    peaks:
        index of peak
    fig:
    labels:
        labels to added to figure
    label_mode:
        if int, it will do one of the following
        if Sequence[int], it will do all listed
        0: None (hover over for plotly)
        1: plot
        2: legend (bounds only)
    mode:
        if int, it will do one of the following
        if Sequence[int], it will do all listed
        0: shaded area (bounds only)
        1: trace of peak signal (bounds only)
        2: range of bounds (bounds only)
        3: max point added (peak index will be used, otherwise max value in bound range will be used)
    normalize:
        0: no normalization
        1: normalize by height
        2: normalize by area
    plot_kwargs

    Returns
    -------

    """
    plotting_lib = global_config.get_plotting_lib() if fig is None else global_config.get_plotting_lib_from_fig(fig)
    if bounds is None and peaks is None:
        raise ValueError("'bounds' and 'peaks' both cannot be 'None'.")
    if bounds is not None and len(bounds.shape) != 2 and bounds.shape[1] != 2:
        raise ValueError("'bound_index' must be 2-dimensional with shape [n,2].")
    if peaks is not None and len(peaks.shape) != 1:
        raise ValueError("'bound_index' must be 2-dimensional with shape [n,2].")

    mode = [mode] if isinstance(mode, int) else mode
    label_mode = [label_mode] if isinstance(label_mode, int) else label_mode

    x = signal_.x
    if normalize == 1:
        y = signal_.y_normalized_by_max()
    elif normalize == 2:
        y = signal_.y_normalized_by_area()
    else:
        y = signal_.y

    labels = "peaks" if labels is None and 0 not in label_mode else labels

    if (3 in mode or 1 in label_mode) and peaks is None:
        peaks = np.array([np.argmax(y[b[0]:b[1]]) + b[0] for b in bounds])
    if labels is None:
        labels = [f"peak: {i}" for i in np.arange(len(bounds))]

    plot_kwargs = copy.copy(plot_kwargs) or {}
    if plotting_lib == global_config.PLOTTING_LIBRARIES.plotly:
        from chem_analysis.plotting.plotly_plots.plotly_peaks import plotly_peaks, plotly_bounds
        if 3 in mode or 1 in label_mode:
            fig = plotly_peaks(fig, plot_kwargs, x, y, peaks, labels, mode, label_mode)
        if bounds is not None:
            fig = plotly_bounds(fig, plot_kwargs, x, y, bounds, labels, mode, label_mode)
        return fig

    if plotting_lib == global_config.PLOTTING_LIBRARIES.matplotlib:
        from chem_analysis.plotting.matplotlib_plots.matplotlib_peaks import matplotlib_peaks, matplotlib_bounds
        if 3 in mode or 1 in label_mode:
            fig = matplotlib_peaks(fig, plot_kwargs, x, y, peaks, labels, mode, label_mode)
        if bounds is not None:
            fig = matplotlib_bounds(fig, plot_kwargs, x, y, bounds, labels, mode, label_mode)
        return fig

    if plotting_lib == global_config.PLOTTING_LIBRARIES.pygraphqt:
        raise NotImplementedError("not done yet")

    raise RuntimeError('Should never happen')


def calibration(
        calibration_,
        *,
        fig=None,
        plot_kwargs: dict | None = None,
):
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_calibration import plotly_calibration
            return plotly_calibration(calibration_, plot_kwargs, fig)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def baseline(
        baseline_: Baseline | Processor,
        *,
        fig=None,
        plot_kwargs: dict | None = None,
):
    plot_kwargs = copy.copy(plot_kwargs) or {}
    if isinstance(baseline_, Processor):
        if not baseline_.processed:
            raise ValueError(
                "Processer has not been used to process a signal yet. Run the 'Processor.run_xy(signal)' method.")
        baselines = [method_ for method_ in baseline_.methods if
                     isinstance(method_, Baseline) and method_.baseline is not None]
        if len(baselines) == 0:
            raise ValueError("No Baseline methods detected. Ensure 'Baseline.save_result' attribute is set to 'True'. ")
        baseline_ = baselines[-1]  # only look at the first one

    if baseline_.baseline is None:
        raise RuntimeError("No baseline detected.\n The 'Baseline.save_result' attribute was likely not set to 'True'.")

    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_baseline import plotly_baseline
            return plotly_baseline(baseline_, plot_kwargs, fig)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()
