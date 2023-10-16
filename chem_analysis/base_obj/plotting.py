import enum
import logging

import numpy as np
import plotly.graph_objs as go  # TODO: move somewhere else

from chem_analysis.config import global_config
from chem_analysis.utils.plot_format import get_plot_color, bold_in_html
from chem_analysis.base_obj.peak import Peak
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.base_obj.signal_array import SignalArray

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
        self.x_range: slice | None = None
        self.y_range: slice | None = None
        self.x_index_range: slice | None = None
        self.y_index_range: slice | None = None

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
                
    def apply_default_formatting(self) -> tuple[dict, dict, dict]:
        layout = {}
        layout_xaxis = {}
        layout_yaxis = {}

        if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
            layout = {
                # "autosize": False,
                # "width": 1200,
                # "height": 600,
                "font": dict(family="Arial", size=18, color="black"),
                "plot_bgcolor": "white",
                # "legend": {"x": 0.05, "y": 0.95}
            }

            layout_xaxis = {
                "tickprefix": "<b>",
                "ticksuffix": "</b>",
                "showline": True,
                "linewidth": 5,
                # "mirror": True,
                "linecolor": 'black',
                "ticks": "outside",
                "tickwidth": 4,
                "showgrid": False,
                "gridwidth": 1,
                "gridcolor": 'lightgray'
            }

            layout_yaxis = {
                "tickprefix": "<b>",
                "ticksuffix": "</b>",
                "showline": True,
                "linewidth": 5,
                "mirror": True,
                "linecolor": 'black',
                "ticks": "outside",
                "tickwidth": 4,
                "showgrid": False,
                "gridwidth": 1,
                "gridcolor": 'lightgray'
            }

        return layout, layout_xaxis, layout_yaxis

    ## signal ########################################################################################################
    def get_x_index_range(self, x: np.ndarray) -> slice:
        if self.x_index_range is not None:
            return self.x_index_range
        if self.x_range is not None:
            start = np.argmin(np.abs(x - self.x_range.start))
            end = np.argmin(np.abs(x - self.x_range.stop))
            if start > end:
                start, end = end, start
            if self.x_range.step is None:
                return slice(start, end)

            avg_step = np.mean(x[1:] - x[:-1])
            step = np.ceil(self.x_range.step / avg_step)
            return slice(start, end, step)

        return slice(0, -1)

    def get_signal_xy(self, signal: Signal) -> tuple[np.ndarray, np.ndarray]:
        y = self.get_y_from_signal(signal)
        slice_ = self.get_x_index_range(signal.x)
        return signal.x[slice_], y[slice_]

    def get_y_index_range(self, y: np.ndarray) -> slice:
        if self.y_index_range is not None:
            return self.y_index_range
        if self.y_range is not None:
            start = np.argmin(y - self.y_range.start)
            end = np.argmin(y - self.y_range.stop)
            if self.y_range.step is None:
                return slice(start, end)
            
            avg_step = np.mean(y[1:] - y[:-1])
            step = np.ceil(self.y_range.step/avg_step)
            return slice(start, end, step)

        return slice(0, -1)

    def get_y(self, signals: list[Signal], y: np.ndarray) -> list[Signal]:
        return signals[self.get_y_index_range(y)]

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

    def get_y_from_signal(self, signal: Signal) -> np.ndarray:
        if self.normalize == NormalizationOptions.AREA:
            return signal.y_normalized_by_area
        if self.normalize == NormalizationOptions.PEAK_HEIGHT:
            return signal.y_normalized_by_peak_max
        return signal.y

    ## peak ########################################################################################################
    @property
    def peak_show(self) -> bool:
        vars_ = [getattr(self, k) for k in vars(self) if k.startswith("peak_show")]
        return any(vars_)

    @peak_show.setter
    def peak_show(self, value: bool):
        for k in vars(self):
            if k.startswith("peak_show"):
                setattr(self, k, value)

    def get_peak_bound_group(self, peak: Peak) -> str:
        ...
    
    def get_peak_bound_color(self, peak: Peak) -> str:
        ...


## Plotting general ###################################################################################################
def plot_signal_array_overlap(
        array: SignalArray,
        *,
        config: PlotConfig = None,
        **kwargs
):
    if array.number_of_signals == 0:
        raise ValueError("No signals to plot.")

    _default_config = kwargs.pop("_default_config")
    if config is None:
        config = _default_config()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)
    config.set_attrs_from_kwargs(**kwargs)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array(array, config)

    raise NotImplementedError()


def plot_signal_array_3D(
        array: SignalArray,
        *,
        config: PlotConfig = None,
        **kwargs
):
    if len(array.signals) == 0:
        raise ValueError("No signals to plot.")

    _default_config = kwargs.pop("_default_config")
    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)
    config.set_attrs_from_kwargs(**kwargs)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array_3D(array, config)

    raise NotImplementedError()


def plot_signal_array_surface(
        array: SignalArray,
        *,
        config: PlotConfig = None,
        **kwargs
):
    if len(array.signals) == 0:
        raise ValueError("No signals to plot.")

    _default_config = kwargs.pop("_default_config")
    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)
    config.set_attrs_from_kwargs(**kwargs)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array_surface(array, config)

    raise NotImplementedError()


def plot_chromatogram(
        chromatogram: Chromatogram,
        *,
        config: PlotConfig = None,
        **kwargs
):
    if len(chromatogram.signals) == 0:
        raise ValueError("No signals to plot.")

    _default_config = kwargs.pop("_default_config")
    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.DIVERSE
    config.set_color_count(chromatogram.number_of_signals)
    config.set_attrs_from_kwargs(**kwargs)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_chromatogram(chromatogram, config)

    raise NotImplementedError()


def plot_signal(
        signal: Signal,
        *,
        config: PlotConfig = None,
        **kwargs
):
    _default_config = kwargs.pop("_default_config")
    if config is None:
        config = PlotConfig()
    config.set_attrs_from_kwargs(**kwargs)
    
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        fig = go.Figure()
        if config.peak_show:
            for peak in signal.peaks:
                plotly_add_peaks(fig, peak, config)
        plotly_signal(fig, signal, config)
        plotly_layout(fig, signal, config)
        return fig

    raise NotImplementedError()


## plotly ##
#######################################################################################################################
#######################################################################################################################
def plotly_signal_array(array: SignalArray, config: PlotConfig):
    fig = go.Figure()

    signals = config.get_y(array.signals, array.y)
    for signal in signals:
        plotly_signal(fig, signal, config)
    plotly_layout(fig, signal, config)
    return fig


def plotly_chromatogram(chromatogram: Chromatogram, config: PlotConfig):
    fig = go.Figure()
    for signal in chromatogram.signals:
        plotly_signal(fig, signal, config)
    plotly_layout(fig, signal, config)
    return fig


def plotly_layout(fig: go.Figure, signal: Signal, config: PlotConfig):
    if config.default_formatting:
        layout, layout_xaxis, layout_yaxis = config.apply_default_formatting()
        fig.update_layout(**layout)
        fig.update_xaxes(**layout_xaxis)
        fig.update_yaxes(**layout_yaxis)

    # label axis
    fig.layout.xaxis.title = bold_in_html(config.get_x_label(signal))
    fig.layout.yaxis.title = bold_in_html(config.get_y_label(signal))


def plotly_signal(fig: go.Figure,
                  signal: Signal,
                  config: PlotConfig
                  ) -> go.Figure:
    x, y = config.get_signal_xy(signal)
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=signal.name,
            connectgaps=config.signal_connect_gaps,
            legendgroup=config.get_signal_group(signal),
            line={"color": config.get_signal_color(signal), "width": config.signal_line_width}
        )
    )

    return fig


def plotly_add_peaks(fig: go.Figure, peak: Peak, config: PlotConfig):
    if config.peak_show_shade:
        plotly_add_peak_shade(fig, peak, config)

    # if config.peak_show_trace:
    #     plotly_add_peak_trace()

    if config.peak_show_bounds:
        plotly_add_peak_bounds(fig, peak, config)

    if config.peak_show_max:
        plotly_add_peak_max(fig, peak, config)


def plotly_add_peak_shade(fig: go.Figure, peak: Peak, config: PlotConfig):
    """ Plots the shaded area for the peak. """
    fig.add_trace(go.Scatter(
        x=peak.x,
        y=peak.y,
        mode="lines",
        fill='tozeroy',
        line={"width": 0},
        showlegend=False,
        legendgroup=config.get_peak_bound_group(peak)
    ))


def plotly_add_peak_max(fig: go.Figure, peak: Peak, config: PlotConfig):
    """ Plots peak name at max. """
    fig.add_trace(go.Scatter(
        x=[peak.max_location],
        y=[peak.max_value],
        mode="text",
        marker={"size": config.peak_marker_size},
        text=[f"{peak.id_}"],
        textposition="top center",
        showlegend=False,
        legendgroup=config.get_peak_bound_group(peak)
    ))


def plotly_add_peak_bounds(fig: go.Figure, peak: Peak, config: PlotConfig):
    """ Adds bounds at the bottom of the plot_add_on for peak area. """ 
    if config.normalize == config.normalizations.AREA:
        bound_height = np.max(peak.parent.y_normalized_by_area) * config.peak_bound_height
    elif config.normalize == config.normalizations.PEAK_HEIGHT:
        bound_height = config.peak_bound_height
    else:
        bound_height = np.max(peak.parent.y_normalized_by_peak_max) * config.peak_bound_height
        
    # bounds
    fig.add_trace(go.Scatter(
        x=[peak.low_bound_location, peak.low_bound_location],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": config.get_peak_bound_color(peak)},
        showlegend=False,
        legendgroup=config.get_peak_bound_group(peak)
    ))
    fig.add_trace(go.Scatter(
        x=[peak.high_bound_location, peak.high_bound_location],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": config.get_peak_bound_color(peak)},
        showlegend=False,
        legendgroup=config.get_peak_bound_group(peak)
    ))
    fig.add_trace(go.Scatter(
        x=[peak.low_bound_location, peak.high_bound_location],
        y=[0, 0],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": config.get_peak_bound_color(peak)},
        showlegend=False,
        legendgroup=config.get_peak_bound_group(peak)
    ))


def layout_3D(fig: go.Figure):
    fig.update_layout(scene=dict(
        yaxis_title='rxn time (sec)',
        xaxis_title='wavenumber (cm-1)',
        zaxis_title='signal'),
    )


def plotly_signal_array_3D(array: SignalArray, config: PlotConfig) -> go.Figure:
    fig = go.Figure()

    for i, t in enumerate(times):
        fig.add_trace(
            go.Scatter3d(
                x=wavenumber,
                y=t * np.ones_like(wavenumber),
                z=data[i, :],
                mode="lines",
                line={"color": "black"},
                legendgroup="lines",
                showlegend=False if i != 0 else True,
            )
        )

    return fig


def plotly_signal_array_surface(array: SignalArray, config: PlotConfig) -> go.Figure:
    fig = go.Figure()

    fig.add_trace(
        go.Surface(
            x=wavenumber,
            y=times,
            z=data,
            legendgroup="surface",
            showlegend=True,
            showscale=False
        )
    )

    return fig

## matplotlib ##
#######################################################################################################################
#######################################################################################################################
