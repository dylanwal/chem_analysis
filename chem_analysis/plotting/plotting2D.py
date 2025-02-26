import logging
import copy
from typing import Sequence

from chem_analysis.config import global_config
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.analysis.peak_result import ResultPeaks2D


logger = logging.getLogger(__name__)


def signal2D_contour(
        signal_: Signal2D,
        *,
        fig=None,
        raw: bool = False,
        plot_kwargs: dict | None = None,
):
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_countour import plotly_contour
            return plotly_contour(signal_, fig, raw, plot_kwargs)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_surface(
        signal_: Signal2D,
        *,
        fig=None,
        raw: bool = False,
        plot_kwargs: dict | None = None
):
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_surface import plotly_surface
            return plotly_surface(signal_, fig, raw, plot_kwargs)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_overlap_signals(
        signal_: Signal2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        raw: bool = False,
        plot_kwargs: dict | None = None,
):
    """

    Parameters
    ----------
    signal_
    slices:
        None: all slices
        int: the index of the slice (y-axis)
        Sequence[int]: the indexes of the slices (y-axis)
        slice:
    fig
    raw

    Returns
    -------

    """
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_slices import plotly_overlap_signals
            return plotly_overlap_signals(signal_, slices, fig, raw, plot_kwargs)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_stack_signals(
        signal_: Signal2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        raw: bool = False,
        plot_kwargs: dict | None = None,
):
    """

    Parameters
    ----------
    signal_
    slices:
        None: all slices
        int: the index of the slice (y-axis)
        Sequence[int]: the indexes of the slices (y-axis)
        slice:
    fig
    raw
    plot_kwargs

    Returns
    -------

    """
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_slices import plotly_stack_signals
            return plotly_stack_signals(signal_, slices, fig, raw, plot_kwargs)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_slices_peaks(
        peaks_: ResultPeaks2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        mode: int | Sequence[int] = 0,
        normalize: int = 0,
        plot_kwargs: dict | None = None,
):
    """

    Parameters
    ----------
    peaks_
    slices
    fig:

    mode:
        if int, it will do one of the following
        if Sequence[int], it will do the following that are listed
        0: shaded area
        1: trace of peak signal
        2: integration bounds
        3: max point added
    normalize:
        0: no normalization
        1: normalize by height
        2: normalize by area
    plot_kwargs:

    Returns
    -------

    """
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_slices import plotly_2D_slices_peaks
            return plotly_2D_slices_peaks(peaks_, slices, fig, mode, normalize, plot_kwargs)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()
