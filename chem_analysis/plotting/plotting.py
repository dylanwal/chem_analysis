import copy
import logging
from typing import Sequence

import numpy as np

from chem_analysis.config import global_config
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.processing.processor import Processor
from chem_analysis.processing.processing_method import Baseline
from chem_analysis.analysis.peak import Peak, PeakData
import chem_analysis.plotting.plotting_primitives as primitives

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
        peaks: Sequence[PeakData],
        *,
        fig=None,
        label_mode: int | Sequence[int] = 0,
        mode: int | Sequence[int] = 0,
        normalize: int = 0,
        colors: str | Sequence[str] | None = None,
):
    """

    Parameters
    ----------
    peaks:
    fig:
    label_mode:
        0: None
        1: legend
    mode:
        if int, it will do one of the following
        if Sequence[int], it will do all listed
        0: shaded area
        1: trace of peak signal
        2: range of bounds
        3: max point added
    normalize:
        0: no normalization
        1: normalize by height
        2: normalize by area
    colors:
        None: default of plotting software (all peaks same color)
        'multi,{color-scheme}': default of plotting software (each peak is a different color)
        single value: all peaks will have same color
        Sequence: must be same length as 'peaks'

    Returns
    -------

    """
    plotting_lib = global_config.get_plotting_lib() if fig is None else global_config.get_plotting_lib_from_fig(fig)

    if not all(peak.parent is peaks[0].parent for peak in peaks):
        raise ValueError("All peaks must have same parent per call.")
    if isinstance(colors, str) and colors.startswith("multi"):
        if ',' in colors:
            style = colors.split(',')[1]
        else:
            style = None

        if plotting_lib == global_config.PLOTTING_LIBRARIES.plotly:
            from chem_analysis.plotting.plotly_plots.plotly_colors import get_n_colors
        else:
            raise NotImplementedError("not done yet")
        colors = get_n_colors(len(peaks), style)

    if isinstance(colors, Sequence) and len(colors) != len(peaks):
        raise ValueError("The length of `colors` and `peaks` must match.")


    mode = [mode] if isinstance(mode, int) else mode
    label_mode = [label_mode] if isinstance(label_mode, int) else label_mode

    # assumes all peaks have the sample parent
    if normalize == 1:
        multiplier = peaks[0].parent.y_normalized_by_max().max() /peaks[0].parent.y.max()
    elif normalize == 2:
        multiplier = peaks[0].parent.y_normalized_by_area().max() /peaks[0].parent.y.max()
    else:
        multiplier = 1

    container = primitives.DrawingContainer()
    for i, peak in enumerate(peaks):
        color = colors[i] if isinstance(colors, Sequence) else colors
        x_fill = peak.x
        y_fill = peak.y * multiplier
        if peak.to_zero:
            # draw the shape down to zero
            x_fill = np.concatenate((x_fill, [x_fill[-1], x_fill[0]]))
            y_fill = np.concatenate((y_fill, [0, 0]))
        else:
            # close the loop
            x_fill = np.concatenate((x_fill, [x_fill[0]]))
            y_fill = np.concatenate((y_fill, [y_fill[0]]))
        if 0 in mode:
            fill = primitives.Fill(x_fill, y_fill, color)
            container.add_objects(fill)
        if 1 in mode:
            line = primitives.Line(x_fill, y_fill, color)
            container.add_objects(line)
        if 2 in mode:
            bound_height = (peak.properties.max_y - peak.properties.min_y) * 0.06
            # side vertical lines
            lines = (
                primitives.Line(
                    np.array([peak.x[0], peak.x[0]]),
                    np.array([peak.y[0]-bound_height / 2, peak.y[0] + bound_height / 2]),
                    color
                ),
                primitives.Line(
                    np.array([peak.x[-1], peak.x[-1]]),
                    np.array([peak.y[-1]-bound_height / 2, peak.y[-1] + bound_height / 2]),
                    color
                ),
                primitives.Line(
                    np.array([peak.x[0], peak.x[-1]]),
                    np.array([peak.y[0], peak.y[-1]]),
                    color
                ),
            )
            container.add_objects(lines)
        if mode == 3:
            marker = primitives.Dot(peak.properties.max_x, peak.properties.max_y, color)
            container.add_objects(marker)

    if plotting_lib == global_config.PLOTTING_LIBRARIES.plotly:
        from chem_analysis.plotting.plotly_plots.plotly_container import draw_container
        fig = draw_container(fig, container, label_mode)
        return fig

    if plotting_lib == global_config.PLOTTING_LIBRARIES.matplotlib:
        raise NotImplementedError("not done yet")

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
    plotting_lib = global_config.get_plotting_lib() if fig is None else global_config.get_plotting_lib_from_fig(fig)

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

    if plotting_lib == global_config.PLOTTING_LIBRARIES.plotly:
            from chem_analysis.plotting.plotly_plots.plotly_baseline import plotly_baseline
            return plotly_baseline(baseline_, plot_kwargs, fig)

    raise NotImplementedError()
