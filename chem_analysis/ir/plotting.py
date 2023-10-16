

from chem_analysis.ir.ir_signal import IRSignal
from chem_analysis.ir.ir_array import IRArray
import chem_analysis.base_obj.plotting as base_plotting


class PlotConfigIR(base_plotting.PlotConfig):
    ...


def plot_signal_array_overlap(
        array: IRArray,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal_array_overlap(array, config=config, _default_config=PlotConfigIR, **kwargs)


def plot_signal_array_3D(
        array: IRArray,
        *,
        config: PlotConfigIR = None,
        **kwargs
):
    return base_plotting.plot_signal_array_3D(array, config=config, _default_config=PlotConfigIR, **kwargs)


def plot_signal_array_surface(
        array: IRArray,
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
