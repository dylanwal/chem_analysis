
from chem_analysis.config import global_config
from chem_analysis.sec.sec_signal import SECSignal
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.base_obj.plotting import PlotConfig
from chem_analysis.analysis.boundary_detection.boundary_detection import ResultPeakBound


def plot_signal(signal: SECSignal, *, fig=None, config=PlotConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        from chem_analysis.sec.plotting_plotly import plotly_signal
        return plotly_signal(fig, signal, config)

    raise NotImplementedError()


def plot_peaks(peaks: ResultPeakBound, *, fig=None, config=PlotConfig()):
    from chem_analysis.base_obj.plotting import plot_peaks as plot_peaks_base
    return plot_peaks_base(peaks, fig=fig, config=config)


def plot_calibration(calibration: SECCalibration, *, fig=None, config=PlotConfig()):
    from chem_analysis.sec.plotting_plotly import plotly_sec_calibration
    return plotly_sec_calibration(calibration, fig=fig, config=config)
