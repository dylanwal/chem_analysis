import logging
from typing import Sequence

from chem_analysis.config import global_config
from chem_analysis.plotting.config import PlotConfig
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.analysis.peak_result import ResultPeaks2D


logger = logging.getLogger(__name__)


def signal2D_contour(
        signal_: Signal2D,
        *,
        fig=None,
        config: PlotConfig | None = None,
        raw: bool = False,
        plot_kwargs: dict | None = None,
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_countour import plotly_contour
            return plotly_contour(signal_, fig, config, raw, plot_kwargs)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_surface(
        signal_: Signal2D,
        *,
        fig=None,
        config: PlotConfig | None = None,
        raw: bool = False
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_surface import plotly_surface
            return plotly_surface(signal_, fig, config, raw)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_slices(
        signal_: Signal2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        config: PlotConfig | None = None,
        raw: bool = False,
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
    config
    raw

    Returns
    -------

    """
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_slices import plotly_slices
            return plotly_slices(signal_, slices, fig, config, raw)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_slices_separated(
        signal_: Signal2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        config: PlotConfig | None = None,
        raw: bool = False,
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
    config
    raw

    Returns
    -------

    """
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_slices import plotly_slices
            return plotly_slices(signal_, slices, fig, config, raw)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def signal2D_slice_peaks(
        peaks_: ResultPeaks2D,
        slices: None | int | Sequence[int] | slice = None,
        *,
        fig=None,
        config: PlotConfig | None = None,
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_2D_peaks import plotly_2D_slice_peaks
            return plotly_2D_slice_peaks(peaks_, slices, fig, config)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


