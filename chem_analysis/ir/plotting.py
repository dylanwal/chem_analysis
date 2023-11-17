
from chem_analysis.config import global_config
from chem_analysis.ir.ir_signal import IRSignal
from chem_analysis.ir.ir_array import IRSignalArray
import chem_analysis.base_obj.plotting as base_plotting
from chem_analysis.base_obj.plotting import PlotConfig
from chem_analysis.analysis.boundary_detection.boundary_detection import ResultPeakBound


class PlotConfigIR(base_plotting.PlotConfig):

    def apply_default_formatting(self) -> tuple[dict, dict, dict]:
        layout, layout_xaxis, layout_yaxis = super().apply_default_formatting()
        layout_xaxis["autorange"] = "reversed"
        return layout, layout_xaxis, layout_yaxis


def plot_signal(signal: IRSignal, *, fig=None, config=PlotConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        from chem_analysis.base_obj.plotting_plotly import plotly_signal
        return plotly_signal(fig, signal, config)

    raise NotImplementedError()


def plot_peaks(peaks: ResultPeakBound, *, fig=None, config=PlotConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        from chem_analysis.base_obj.plotting_plotly import plotly_peaks
        return plotly_peaks(fig, peaks, config)

    raise NotImplementedError()
