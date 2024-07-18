from typing import Sequence
import copy
import logging

from chem_analysis.config import global_config
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.processing.processing_method import Baseline
from chem_analysis.analysis.peak_result import ResultPeaks

logger = logging.getLogger(__name__)


def signal(
        signal_: Signal,
        *,
        fig=None,
        raw: bool = False,
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
    raw:
        True: x_
    normalize:
        0: no normalization
        1: normalize by height
        2: normalize by area
    plot_kwargs:
        arg passed directly to plotting function

    Returns
    -------

    """
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_signal import plotly_signal
            return plotly_signal(signal_, plot_kwargs, fig, raw, normalize)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def peaks(
        peaks_: ResultPeaks,
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
    fig
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
    plot_kwargs

    Returns
    -------

    """
    plot_kwargs = copy.copy(plot_kwargs) or {}
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_peaks import plotly_peaks
            return plotly_peaks(peaks_, plot_kwargs, fig, mode, normalize)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def calibration(
        calibration_: SECCalibration,
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
        baseline_: Baseline | Signal,
        *,
        fig=None,
        plot_kwargs: dict | None = None,
):
    plot_kwargs = copy.copy(plot_kwargs) or {}
    if isinstance(baseline_, Signal):
        baselines = [method_ for method_ in baseline_.processor.methods if isinstance(method_, Baseline)]
        if len(baselines) == 0:
            raise ValueError("No Baseline methods detected.")
        baseline_ = baselines[0]

    if baseline_.baseline is None:
        raise RuntimeError("No baseline detected.\nEither the processing method has not been run yet (call Signal.x to force processing) or"
                         "the 'Baseline.save_result' attribute was not set to 'True'.")

    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_baseline import plotly_baseline
            return plotly_baseline(baseline_, plot_kwargs, fig)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()
