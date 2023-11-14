
from chem_analysis.ir.ir_signal import IRSignal
from chem_analysis.ir.ir_array import IRSignalArray
import chem_analysis.base_obj.plotting as base_plotting


class PlotConfigIR(base_plotting.PlotConfig):

    def apply_default_formatting(self) -> tuple[dict, dict, dict]:
        layout, layout_xaxis, layout_yaxis = super().apply_default_formatting()
        layout_xaxis["autorange"] = "reversed"
        return layout, layout_xaxis, layout_yaxis


def plot_signal_array_overlap(
        array: IRSignalArray,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal_array_overlap(array, config=config, _default_config=PlotConfigIR, **kwargs)


def plot_signal_array_3D(
        array: IRSignalArray,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal_array_3D(array, config=config, _default_config=PlotConfigIR, **kwargs)


def plot_signal_array_surface(
        array: IRSignalArray,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal_array_surface(array, config=config, _default_config=PlotConfigIR, **kwargs)


def plot_signal(
        signal: IRSignal,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal(signal, config=config, _default_config=PlotConfigIR, **kwargs)
