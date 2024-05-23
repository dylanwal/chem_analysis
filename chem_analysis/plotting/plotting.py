import logging

from chem_analysis.config import global_config
from chem_analysis.plotting.config import PlotConfig
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_discrete import SignalDiscrete
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.processing.processing_method import Baseline
from chem_analysis.analysis.peak_result import ResultPeaks

logger = logging.getLogger(__name__)


def signal(
        signal_: Signal | SignalDiscrete,
        *,
        fig=None,
        config: PlotConfig | None = None,
        raw: bool = False
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_signal import plotly_signal
            return plotly_signal(signal_, fig, config, raw)
        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def peaks(
        peaks_: ResultPeaks,
        *,
        fig=None,
        config: PlotConfig | None = None,
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_peaks import plotly_peaks
            return plotly_peaks(peaks_, fig, config)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def calibration(
        calibration_: SECCalibration,
        *,
        fig=None,
        config: PlotConfig | None = None,
):
    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_calibration import plotly_calibration
            return plotly_calibration(calibration_, fig, config)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()


def get_baseline_from_signal(signal_: Signal):
    baselines = [method_ for method_ in signal_.processor.methods if isinstance(method_, Baseline)]
    if len(baselines) == 0:
        raise ValueError("No Baseline methods detected.")
    return baselines[0]


def baseline(
        baseline_: Baseline | Signal,
        *,
        fig=None,
        config: PlotConfig | None = None,
):
    if isinstance(baseline_, Signal):
        baseline_ = get_baseline_from_signal(baseline_)
    if baseline_.baseline is None:
        raise RuntimeError("No baseline detected.\nEither the processing method has not been run yet (call Signal.x to force processing) or"
                         "the 'Baseline.save_result' attribute was not set to 'True'.")

    for option in global_config.get_plotting_options():
        if option == global_config.PLOTTING_LIBRARIES.PLOTLY:
            from chem_analysis.plotting.plotly_plots.plotly_baseline import plotly_baseline
            return plotly_baseline(baseline_, fig, config)

        if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
            pass

        if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
            pass

    raise NotImplementedError()

#
# def signal2d_dynamic(
#         array_: Signal2D,
#         *,
#         config=PlotConfig()
# ):
#     config = config or PlotConfig()
#     for option in global_config.get_plotting_options():
#         if isinstance(calibration, SECCalibration):
#             pass
#
#         if option == global_config.PLOTTING_LIBRARIES.MATPLOTLIB:
#             pass
#
#         if option == global_config.PLOTTING_LIBRARIES.PYGRAPHQT:
#             from chem_analysis.plotting.qt_plots.qt_array import qt_array
#             return qt_array(array_)
#
#     raise NotImplementedError()

# def signal2d_overlap(
#         array: SignalTimeSeries,
#         *,
#         config: PlotConfig = None,
#         **kwargs
# ):
#     if array.number_of_signals == 0:
#         raise ValueError("No signals to plot.")
#
#     _default_config = kwargs.pop("_default_config")
#     if config is None:
#         config = _default_config()
#         config.signal_color = SignalColorOptions.SINGLE
#     config.set_color_count(array.number_of_signals)
#     config.set_attrs_from_kwargs(**kwargs)
#
#     if global_config.plotting_library == global_config.PLOTTING_LIBRARIES.PLOTLY:
#         return plotly_signal_array(array, config)
#
#     raise NotImplementedError()
#
#
# def plot_signal_array_3D(
#         array: SignalTimeSeries,
#         *,
#         config: PlotConfig = None,
#         **kwargs
# ):
#     if len(array.signals) == 0:
#         raise ValueError("No signals to plot.")
#
#     _default_config = kwargs.pop("_default_config")
#     if config is None:
#         config = _default_config()
#         config.signal_color = SignalColorOptions.SINGLE
#     config.set_color_count(array.number_of_signals)
#     config.set_attrs_from_kwargs(**kwargs)
#
#     if global_config.plotting_library == global_config.PLOTTING_LIBRARIES.PLOTLY:
#         return plotly_signal_array_3D(array, config)
#
#     raise NotImplementedError()
#
#
# def plot_signal_array_surface(
#         array: SignalTimeSeries,
#         *,
#         config: PlotConfig = None,
#         **kwargs
# ):
#     if len(array.signals) == 0:
#         raise ValueError("No signals to plot.")
#
#     _default_config = kwargs.pop("_default_config")
#     if config is None:
#         config = _default_config()
#         config.signal_color = SignalColorOptions.SINGLE
#     config.set_color_count(array.number_of_signals)
#     config.set_attrs_from_kwargs(**kwargs)
#
#     if global_config.plotting_library == global_config.PLOTTING_LIBRARIES.PLOTLY:
#         return plotly_signal_array_surface(array, config)
#
#     raise NotImplementedError()
#
#
# def plot_chromatogram(
#         chromatogram: Chromatogram,
#         *,
#         config: PlotConfig = None,
#         **kwargs
# ):
#     if len(chromatogram.signals) == 0:
#         raise ValueError("No signals to plot.")
#
#     _default_config = kwargs.pop("_default_config")
#     if config is None:
#         config = _default_config()
#         config.signal_color = SignalColorOptions.DIVERSE
#     config.set_color_count(chromatogram.number_of_signals)
#     config.set_attrs_from_kwargs(**kwargs)
#
#     if global_config.plotting_library == global_config.PLOTTING_LIBRARIES.PLOTLY:
#         return plotly_chromatogram(chromatogram, config)
#
#     raise NotImplementedError()
