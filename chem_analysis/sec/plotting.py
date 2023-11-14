
import numpy as np

from chem_analysis.config import global_config
from chem_analysis.sec.sec_signal import SECSignal
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.base_obj.plotting import PlotConfig
from chem_analysis.analysis.boundry_detection.bound_detection import ResultPeakBound


def plot_sec_signal(signal: SECSignal, *, fig=None, config=PlotConfig()):
    from chem_analysis.base_obj.plotting import plot_signal
    return plot_signal(signal, fig=fig, config=config)


def plot_sec_peaks(peaks: ResultPeakBound, *, fig=None, config=PlotConfig()):
    from chem_analysis.base_obj.plotting import plot_peaks
    return plot_peaks(peaks, fig=fig, config=config)


def plot_sec_calibration(calibration: SECCalibration, *, fig=None, config=PlotConfig()):
    from chem_analysis.sec.plotting_plotly import plotly_sec_calibration
    return plotly_sec_calibration(calibration, fig=fig, config=config)
