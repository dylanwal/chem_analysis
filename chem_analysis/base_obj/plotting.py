import enum
import logging

import numpy as np

from chem_analysis.config import global_config
from chem_analysis.utils.plot_format import get_plot_color, bold_in_html
from chem_analysis.analysis.peak import Peak
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.base_obj.signal_array import SignalArray
from chem_analysis.analysis.boundry_detection.bound_detection import ResultPeakBound


logger = logging.getLogger("plotting")


class NormalizationOptions(enum.Enum):
    NONE = 0
    AREA = 1
    PEAK_HEIGHT = 2


class SignalColorOptions(enum.Enum):
    SINGLE = 0
    DIVERSE = 1


class PlotConfig:
    normalizations = NormalizationOptions
    signal_colors = SignalColorOptions
    
    def __init__(self):
        self.default_formatting: bool = True
        self.title = None
        self.y_label = None
        self.x_label = None
        self.z_label = None

        # signal
        self.normalize: NormalizationOptions = NormalizationOptions.NONE
        # self.signal_color: Callable | str = 'rgb(10,36,204)'
        self.signal_line_width = 2
        self._signal_color_counter = 0
        self._set_color_counter = None
        self._set_color_array_ = None

        # peak
        # self.peak_show_trace: bool = False
        self.peak_show_shade: bool = True
        self.peak_show_bounds: bool = False
        self.peak_show_max: bool = True
        self.signal_connect_gaps: bool = False
        self.signal_color: str | SignalColorOptions | None = None
        self.peak_bound_color: str | None = None
        self.peak_bound_line_width: float = 3
        self.peak_bound_height: float = 0.05  # % of max
        self.peak_marker_size: float = 3
        
    def set_attrs_from_kwargs(self, **kwargs):
        if not kwargs:
            return

        for k, v in kwargs.items():
            if not hasattr(self, k):
                logger.warning(f"'{k} is an invalid plot configuration.")
            setattr(self, k, v)

    ## layout ########################################################################################################
    def get_x_label(self, signal: Signal) -> str:
        if self.x_label is not None:
            return self.x_label
        return signal.x_label

    def get_y_label(self, signal: Signal) -> str:
        if self.y_label is not None:
            return self.y_label
        return signal.y_label

    def get_z_label(self, signal: SignalArray) -> str:
        if self.z_label is not None:
            return self.z_label
        return signal.z_label

    ## signal ########################################################################################################
    def get_signal_color(self, signal: Signal) -> str:
        if self._set_color_counter is None:
            return 'rgb(10,36,204)'

        return self._set_color_array[signal.id_]

    @property
    def _set_color_array(self):
        if self._set_color_array_ is None:
            if self.signal_color is None or self.signal_color is SignalColorOptions.SINGLE:
                self._set_color_array_ = get_plot_color(self._set_color_counter)

        return self._set_color_array_

    def set_color_count(self, number_colors: int):
        self._set_color_counter = number_colors

    def get_signal_group(self, signal: Signal) -> str:
        return f"sig_{signal.id_}"


## Plotting general ###################################################################################################
# def plot_signal_array_overlap(
#         array: SignalArray,
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
#     if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
#         return plotly_signal_array(array, config)
#
#     raise NotImplementedError()
#
#
# def plot_signal_array_3D(
#         array: SignalArray,
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
#     if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
#         return plotly_signal_array_3D(array, config)
#
#     raise NotImplementedError()
#
#
# def plot_signal_array_surface(
#         array: SignalArray,
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
#     if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
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
#     if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
#         return plotly_chromatogram(chromatogram, config)
#
#     raise NotImplementedError()
#

def plot_signal(
        signal: Signal,
        *,
        fig=None,
        config: PlotConfig = None,
        **kwargs
):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        from chem_analysis.base_obj.plotting_plotly import plotly_signal
        return plotly_signal(fig, signal, config)

    raise NotImplementedError()


def plot_peaks(
        peaks: ResultPeakBound,
        *,
        fig=None,
        config: PlotConfig = None,
        **kwargs
):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        from chem_analysis.base_obj.plotting_plotly import plotly_peaks
        return plotly_peaks(fig, peaks, config)

    raise NotImplementedError()
