import copy
from typing import Sequence
import logging

import numpy as np
import matplotlib.pyplot as plt

from chem_analysis.plotting.matplotlib_plots.matplotlib_utils import input_check
from chem_analysis.analysis.peak import PeakContinuous

logger = logging.getLogger(__name__)


def matplotlib_peaks(
        fig: plt.Figure | None,
        plot_kwargs: dict,
        x: np.ndarray,
        y: np.ndarray,
        peak_index: np.ndarray,
        labels: Sequence[str] | None,
        mode: Sequence[int],
        mode_labels: Sequence[int],
) -> plt.Figure:
    plot_kwargs = copy.deepcopy(plot_kwargs)
    fig = input_check(fig)
    ax = fig.axes[0] if fig.axes else fig.add_subplot(111)

    if len(peak_index) == 0:
        logger.warning("No peaks supplied.")
        return fig

    # Set default marker properties
    marker = plot_kwargs.get("marker", {})
    marker_color = marker.get("color", "black")
    marker_style = marker.get("style", "o")  # Default to circle

    # Scatter plot for peak markers
    if 3 in mode:
        ax.scatter(x[peak_index], y[peak_index], color=marker_color, marker=marker_style, label=plot_kwargs.get("name", "peaks"))

    # Add text labels if provided
    if labels is not None and 1 in mode_labels:
        y_offset = (np.max(y[peak_index]) - np.min(y[peak_index])) * 0.03
        for idx, label in zip(peak_index, labels):
            ax.text(x[idx], y[idx] + y_offset, label, ha='center', va='bottom', fontsize=10)

    return fig


def matplotlib_bounds(
        fig: plt.Figure | None,
        plot_kwargs: dict,
        x: np.ndarray,
        y: np.ndarray,
        bounds: np.ndarray | None,
        labels: Sequence[str],
        mode_: Sequence[int],
        label_mode: Sequence[int],
) -> plt.Figure:
    plot_kwargs = copy.deepcopy(plot_kwargs)
    fig = input_check(fig)
    ax = fig.axes[0] if fig.axes else fig.add_subplot(111)

    if len(bounds) == 0:
        logger.warning("No bounds supplied.")
        return fig

    for i, bound in enumerate(bounds):
        if 0 in mode_ or 1 in mode_:
            matplotlib_add_peak_trace(
                ax,
                x[bound[0]:bound[1]],
                y[bound[0]:bound[1]],
                labels[i], mode_, label_mode, plot_kwargs, showlegend=True if i == 0 else False
            )
        if 2 in mode_:
            peak_height = np.max(y)
            matplotlib_add_peak_bounds(ax, x[bound[0]], x[bound[1]], peak_height, plot_kwargs)

    min_b = np.min(bounds)
    max_b = np.max(bounds)
    min_ = np.min(y[min_b:max_b])
    max_ = np.max(list([np.max(y[b[0]:b[1]]) for b in bounds]))
    span = max_ - min_
    ax.set_ylim(min_ - span * 0.05, max_ + span * 0.05)

    if 2 in label_mode:
        ax.legend()

    return fig


def matplotlib_add_peak_trace(
        ax: plt.axes,
        x: np.ndarray,
        y: np.ndarray,
        label: str,
        mode: Sequence[int],
        label_mode: Sequence[int],
        plot_kwargs: dict,
        showlegend: bool = False,
):
    """ Plots the shaded area for the peak. """
    kwargs = {
        "label": label if 2 in label_mode else "bounds",
    }
    plot_kwargs2 = {**plot_kwargs, **kwargs}

    if 0 in mode:
        ax.fill_between(x, y, alpha=0.3, **plot_kwargs2)  # Shaded area
    if 1 in mode:
        ax.plot(x, y, linewidth=plot_kwargs2.get("linewidth", 3), **plot_kwargs2)  # Line plot


def matplotlib_add_peak_bounds(
        ax: plt.axes,
        min_x: int | float,
        max_x: int | float,
        peak_height: int | float,
        plot_kwargs: dict
):
    """ Adds bounds at the bottom of the plot_add_on for peak area. """
    line_kwargs = {"linewidth": 1, "color": "black"} | plot_kwargs
    bound_height = peak_height * 0.06

    # Side vertical lines
    ax.plot([min_x, min_x], [-bound_height / 2, bound_height / 2], **line_kwargs)
    ax.plot([max_x, max_x], [-bound_height / 2, bound_height / 2], **line_kwargs)

    # Horizontal line
    ax.plot([min_x, max_x], [0, 0], **line_kwargs)
