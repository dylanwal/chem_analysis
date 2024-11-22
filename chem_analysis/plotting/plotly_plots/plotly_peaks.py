from typing import Sequence
from logging import getLogger

import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.analysis.peak import PeakContinuous
from chem_analysis.analysis.peak_result import ResultPeaks
from chem_analysis.analysis.ms_analysis.result_search import PeakCompound

logger = getLogger(__name__)


def plotly_peaks(
        peaks: ResultPeaks,
        plot_kwargs: dict,
        fig: go.Figure | None,
        mode: int | Sequence[int] = 0,
        normalize: int = 0,
) -> go.Figure:
    if isinstance(mode, int):
        mode = [mode]
    fig = input_check(fig)
    if len(peaks.peaks) == 0:
        logger.warning('No peaks to added to figure.')
        return fig

    if len(peaks.peaks) > 0 and (not hasattr(peaks.peaks[0], "parent") or not hasattr(peaks.peaks[0], "max_x")):
        raise ValueError(f"Not supported peak type.\n\tpeak type received: {type(peaks[0])}")

    if normalize == 1:
        y_max = np.max(peaks.peaks[0].parent.y_normalized_by_max())
    if normalize == 2:
        y_max = np.max(peaks.peaks[0].parent.y_normalized_by_area())

    plot_kwargs["showlegend"] = plot_kwargs.get("showlegend", True)  # show in legend first peak only
    plot_kwargs["legendgroup"] = plot_kwargs.get("legendgroup", f"peaks ({peaks.peaks[0].parent.name})")
    plot_kwargs["name"] = plot_kwargs.get("name", plot_kwargs["legendgroup"])
    for peak in peaks.peaks:
        y = peak.y
        if normalize > 0:
            y = y / y_max

        if 0 in mode or 1 in mode:
            plotly_add_peak_trace(fig, peak.x, y, get_hover_stats(peak), mode, plot_kwargs)
        if 2 in mode:
            peak_height = np.max(y)
            plotly_add_peak_bounds(fig, peak.min_x, peak.max_x, peak_height, plot_kwargs)
        if 3 in mode:
            peak_height = np.max(y)
            if isinstance(peak, PeakCompound) and peak.compound is not None:
                text = peak.compound.label
            else:
                text = f"peak {peak.id_}"
            plotly_add_peak_max(fig, peak.max_x, peak_height, text, plot_kwargs)

        plot_kwargs["showlegend"] = False

    return fig


def get_hover_stats(peak: PeakContinuous):
    text = []
    if isinstance(peak, PeakCompound) and peak.compound is not None:
        text.append(f"compound: {peak.compound.label}")
    else:
        text.append(f"label: {peak.id_}")
    text += [
        f"span: [{peak.bounds[0]:.2f}, {peak.bounds[1]:.2f}]",
        f"max: {peak.max_y:,.2f} at {peak.max_x:.2f}",
        f"area: {peak.area():,.2f}"
    ]
    return "<br>".join(text)


def plotly_add_peak_trace(
        fig: go.Figure,
        x: np.ndarray,
        y: np.ndarray,
        hover_stats: str,
        mode: Sequence[int],
        plot_kwargs: dict
):
    """ Plots the shaded area for the peak. """
    kwargs = dict(
        x=x,
        y=y,
        mode="lines",
        hovertemplate='<b>%{customdata}</b>',
        customdata=[hover_stats] * len(x),
    )
    plot_kwargs2 = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    if 0 in mode:
        plot_kwargs2["fill"] = plot_kwargs2.get("fill", 'tozeroy')
        plot_kwargs2["line"] = {"width": 0} | plot_kwargs2.get("line", dict())
    if 2 in mode:
        width = plot_kwargs2.get("line", dict()).get("width", 2)
        plot_kwargs2["line"]["width"] = width if width != 0 else 2

    fig.add_scatter(**plot_kwargs2)


def plotly_add_peak_max(fig: go.Figure, max_x: int | float, max_y: int | float, text: str, plot_kwargs: dict):
    """ Plots peak name at max. """
    kwargs = dict(
        x=[max_x],
        y=[max_y],
        mode="text",
        marker={"size": 3},
        text=[text],
        textposition="top center",
    )
    fig.add_scatter(**(kwargs | plot_kwargs))


def plotly_add_peak_bounds(
        fig: go.Figure,
        min_x: int | float,
        max_x: int | float,
        peak_height: int | float,
        plot_kwargs: dict
):
    """ Adds bounds at the bottom of the plot_add_on for peak area. """
    line = {"width": 1, "color": 'rgb(0,0,0)'},
    bound_height = peak_height * 0.06

    # side vertical lines
    fig.add_scatter(
        x=[min_x, max_x],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line=line,
        **plot_kwargs
    )
    fig.add_scatter(
        x=[min_x, max_x],
        y=[-bound_height / 2, bound_height / 2],
        mode="lines",
        line=line,
        **plot_kwargs
    )
    # horizontal line
    fig.add_scatter(
        x=[min_x, max_x],
        y=[0, 0],
        mode="lines",
        line=line,
        **plot_kwargs
    )
