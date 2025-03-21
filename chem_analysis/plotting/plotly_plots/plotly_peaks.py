import copy
from typing import Sequence
from logging import getLogger

import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check

logger = getLogger(__name__)


def plotly_peaks(
        fig: go.Figure | None,
        plot_kwargs: dict,
        x: np.ndarray,
        y: np.ndarray,
        peak_index: np.ndarray,
        labels: Sequence[str] | None,
        mode: Sequence[int],
        label_mode: Sequence[int],
) -> go.Figure:
    plot_kwargs = copy.deepcopy(plot_kwargs)
    fig = input_check(fig)
    if len(peak_index) == 0:
        logger.warning("No peaks supplied.")
        return fig

    plot_kwargs["showlegend"] = plot_kwargs.get("showlegend", True)  # show in legend first peak only
    plot_kwargs["legendgroup"] = plot_kwargs.get("legendgroup", "peaks")
    plot_kwargs["name"] = plot_kwargs.get("name", plot_kwargs["legendgroup"])
    marker = plot_kwargs.get("marker", dict())
    marker['color'] = marker.get("color", "black" if 3 in mode else "rgba(0,0,0,0)")
    plot_kwargs["marker"] = marker
    if 0 in label_mode:
        plot_kwargs["hovertemplate"] = '<b>%{text}</b>'
        plot_kwargs["text"] = labels or [f"peak: {i}" for i in np.arange(len(peak_index))]
    if 1 in label_mode:
        plot_kwargs["text"] = labels
        plot_kwargs["mode"] = "markers + text"
        plot_kwargs["textposition"] = "top center"
    else:
        plot_kwargs["mode"] = "markers"

    fig.add_scatter(x=x[peak_index], y=y[peak_index], **plot_kwargs)
    return fig


def plotly_bounds(
        fig: go.Figure | None,
        plot_kwargs: dict,
        x: np.ndarray,
        y: np.ndarray,
        bounds: np.ndarray | None,
        labels: Sequence[str],
        mode_: Sequence[int],
        label_mode: Sequence[int],
) -> go.Figure:
    plot_kwargs = copy.deepcopy(plot_kwargs)
    fig = input_check(fig)
    if len(bounds) == 0:
        logger.warning("No bounds supplied.")
        return fig

    for i, bound in enumerate(bounds):
        if 0 in mode_ or 1 in mode_:
            plotly_add_peak_trace(
                fig,
                x[bound[0]:bound[1]],
                y[bound[0]:bound[1]],
                labels[i], mode_, label_mode, plot_kwargs, showlegend=True if i == 0 else False)
        if 2 in mode_:
            peak_height = np.max(y)
            plotly_add_peak_bounds(fig, x[bound[0]], x[bound[1]], peak_height, plot_kwargs)

    min_b = np.min(bounds)
    max_b = np.max(bounds)
    min_ = np.min(y[min_b:max_b])
    max_ = np.max(list([np.max(y[b[0]:b[1]]) for b in bounds]))
    span = max_ - min_
    fig.layout.yaxis.range = (min_-span*0.05, max_+span*0.05)

    return fig


def plotly_add_peak_trace(
        fig: go.Figure,
        x: np.ndarray,
        y: np.ndarray,
        label: str,
        mode: Sequence[int],
        label_mode: Sequence[int],
        plot_kwargs: dict,
        showlegend: bool = False,
):
    """ Plots the shaded area for the peak. """
    kwargs = dict(
        x=x,
        y=y,
        mode="lines",
        hovertemplate='<b>[%{x},%{y}]' + f'<br>{label}</b>',
        legendgroup='peaks',
        showlegend=showlegend,
    )
    if 2 in label_mode:
        kwargs["name"] = label
        kwargs["showlegend"] = True
    else:
        kwargs["name"] = "bounds"

    plot_kwargs2 = kwargs | plot_kwargs  # plot_kwargs overwrite kwargs
    if 0 in mode:
        plot_kwargs2["fill"] = plot_kwargs2.get("fill", 'tozeroy')
        plot_kwargs2["line"] = {"width": 0} | plot_kwargs2.get("line", dict())
    if 1 in mode:
        line = plot_kwargs2.get("line", dict())
        line['width'] = line.get("width", 3)
        plot_kwargs2["line"] = line

    fig.add_scatter(**plot_kwargs2)


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
