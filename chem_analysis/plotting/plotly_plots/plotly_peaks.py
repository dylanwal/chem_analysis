
import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.analysis.peak import PeakBounded
from chem_analysis.analysis.peak_result import ResultPeaks
from chem_analysis.analysis.ms_analysis.search_result import PeakCompound


def plotly_peaks(peaks: ResultPeaks, fig: go.Figure | None, config: PlotlyConfig | None) -> go.Figure:
    fig, config = PlotlyConfig.input_check(fig, config)

    # if not isinstance(peaks.peaks[0], PeakBounded):
    #     raise ValueError(f"Not supported peak type.\n\tpeak type received: {type(peaks[0])}")

    for peak in peaks.peaks:
        if isinstance(peak, PeakCompound) and peak.compound is not None:
            label = peak.compound.label
        else:
            label = f"peak {peak.id_}"

        if config.peak_show_shade:
            plotly_add_peak_shade(fig, peak, config, label)
        # if config.peak_show_trace:
        #     plotly_add_peak_trace()
        if config.peak_show_bounds:
            plotly_add_peak_bounds(fig, peak, config, label)
        if config.peak_show_max:
            plotly_add_peak_max(fig, peak, config, label)

    return fig


def plotly_add_peak_shade(fig: go.Figure, peak: PeakBounded, config: PlotlyConfig, label: str):
    """ Plots the shaded area for the peak. """
    fig.add_trace(go.Scatter(
        x=peak.x,
        y=peak.y,
        mode="lines",
        fill='tozeroy',
        line={"width": 0},
        showlegend=True,
        legendgroup=label,
        hovertemplate='<b>%{customdata}</b>',
        customdata=[get_hover_stats(peak)]*len(peak.x),
        name=label
    ))


def get_hover_stats(peak: PeakBounded):
    text = [
        f"id: {peak.id_}",
        f"span: [{peak.low_bound_x:.2f}, {peak.high_bound_x:.2f}]",
        f"max: {peak.max_y:,.2f} at {peak.max_x:.2f}",
        f"area: {peak.area():,.2f}"
    ]
    if isinstance(peak, PeakCompound) and peak.compound is not None:
        text.append(f"compound: {peak.compound.label}")

    return "<br>".join(text)


def plotly_add_peak_max(fig: go.Figure, peak: PeakBounded, config: PlotlyConfig, label: str):
    """ Plots peak name at max. """
    fig.add_trace(go.Scatter(
        x=[peak.max_x],
        y=[peak.max_y],
        mode="text",
        marker={"size": config.peak_marker_size},
        text=[f"{peak.id_}"],
        textposition="top center",
        showlegend=False,
        legendgroup=label
    ))


def plotly_add_peak_bounds(fig: go.Figure, peak: PeakBounded, config: PlotlyConfig, label: str):
    """ Adds bounds at the bottom of the plot_add_on for peak area. """
    if config.normalize == config.NORMALIZATION_OPTIONS.AREA:
        bound_height = np.max(peak.max_y) * config.peak_bound_height
    elif config.normalize == config.NORMALIZATION_OPTIONS.PEAK_HEIGHT:
        bound_height = config.peak_bound_height
    else:
        bound_height = np.max(peak.max_y) * config.peak_bound_height

    # bounds
    fig.add_trace(go.Scatter(
        x=[peak.min_x, peak.min_x],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": 'rgb(0,0,0)'},
        showlegend=False,
        legendgroup=label
    ))
    fig.add_trace(go.Scatter(
        x=[peak.max_x, peak.max_x],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": 'rgb(0,0,0)'},
        showlegend=False,
        legendgroup=label
    ))
    fig.add_trace(go.Scatter(
        x=[peak.min_x, peak.max_x],
        y=[0, 0],
        mode="lines",
        line={"width": config.peak_bound_line_width, "color": 'rgb(0,0,0)'},
        showlegend=False,
        legendgroup=label
    ))
