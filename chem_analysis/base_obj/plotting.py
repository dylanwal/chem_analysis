import enum

import numpy as np
import plotly.graph_objs as go  #TODO: move somewhere else

from chem_analysis.config import global_config
from chem_analysis.utils.plot_format import get_plot_color, bold_in_html
from chem_analysis.base_obj.peak import Peak
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.base_obj.signal_array import SignalArray


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

        # signal
        self.normalize: NormalizationOptions = NormalizationOptions.NONE
        # self.signal_color: Callable | str = 'rgb(10,36,204)'
        self.signal_line_width = 2
        self._signal_color_counter = 0

        # peak
        # self.peak_show_trace: bool = False
        self.peak_show_shade: bool = True
        self.peak_show_bounds: bool = False
        self.peak_show_max: bool = True
        self.signal_connect_gaps: bool = False
        self.signal_color: str | SignalColorOptions = None
        self.peak_bound_color: str | None = None
        self.peak_bound_line_width: float = 3
        self.peak_bound_height: float = 0.05  # % of max
        self.peak_marker_size: float = 3

        self._set_color_count = None

    @property
    def peak_show(self) -> bool:
        vars_ = [getattr(self, k) for k in vars(self) if k.startswith("peak_show")]
        return any(vars_)

    @peak_show.setter
    def peak_show(self, value: bool):
        for k in vars(self):
            if k.startswith("peak_show"):
                setattr(self, k, value)
                
    def apply_default_formatting(self, fig):
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
                "mirror": True,
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
            fig.update_layout(**layout)
            fig.update_xaxes(**layout_xaxis)
            fig.update_yaxes(**layout_yaxis)

    def get_signal_color(self, signal: Signal) -> str:
        if self._set_color_count is None:
            return 'rgb(10,36,204)'

    def get_signal_group(self, signal: Signal) -> str:
        ...

    def set_color_count(self, number_colors: int):
        self._set_color_count = number_colors
        
    def get_peak_bound_group(self, peak: Peak) -> str:
        ...
    
    def get_peak_bound_color(self, peak: Peak) -> str:
        ...


def plot_signal_array_overlap(array: SignalArray, config: PlotConfig = None):
    if len(array.signals) == 0:
        raise ValueError("No signals to plot.")

    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array(array, config)

    raise NotImplementedError()


def plot_signal_array_3D(array: SignalArray, config: PlotConfig = None):
    if len(array.signals) == 0:
        raise ValueError("No signals to plot.")

    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array_3D(array, config)

    raise NotImplementedError()


def plot_signal_array_surface(array: SignalArray, config: PlotConfig = None):
    if len(array.signals) == 0:
        raise ValueError("No signals to plot.")

    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.SINGLE
    config.set_color_count(array.number_of_signals)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_signal_array_surface(array, config)

    raise NotImplementedError()


def plot_chromatogram(chromatogram: Chromatogram, config: PlotConfig = None):
    if len(chromatogram.signals) == 0:
        raise ValueError("No signals to plot.")

    if config is None:
        config = PlotConfig()
        config.signal_color = SignalColorOptions.DIVERSE
    config.set_color_count(chromatogram.number_of_signals)

    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        return plotly_chromatogram(chromatogram, config)

    raise NotImplementedError()


def plot_signal(signal: Signal, config: PlotConfig = PlotConfig()):
    if global_config.plotting_library == global_config.plotting_libraries.PLOTLY:
        fig = go.Figure()
        plotly_signal(fig, signal, config)
        plotly_layout(fig, signal, config)
        return fig

    raise NotImplementedError()


## plotly ##
#######################################################################################################################
#######################################################################################################################
def plotly_signal_array(array: SignalArray, config: PlotConfig):
    fig = go.Figure()
    for signal in array.signals:
        plotly_signal(fig, signal, config)
    plotly_layout(fig, signal, config)


def plotly_chromatogram(chromatogram: Chromatogram, config: PlotConfig):
    fig = go.Figure()
    for signal in chromatogram.signals:
        plotly_signal(fig, signal, config)
    plotly_layout(fig, signal, config)


def plotly_layout(fig: go.Figure, signal: Signal, config: PlotConfig):
    if config.default_formatting:
        config.apply_default_formatting(fig)

    # label axis
    if config.x_label is not None:
        x_label = config.x_label
    else:
        x_label = signal.x_label
    fig.layout.xaxis.title = bold_in_html(x_label)

    if config.y_label is not None:
        y_label = config.y_label
    else:
        y_label = signal.y_label
    fig.layout.yaxis.title = bold_in_html(y_label)


def plotly_signal(fig: go.Figure, signal: Signal, config: PlotConfig) -> go.Figure:
    if config.peak_show:
        for peak in signal.peaks:
            plotly_add_peaks(fig, peak, config)

    fig.add_trace(
        go.Scatter(
            x=signal.x,
            y=signal.y,
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
