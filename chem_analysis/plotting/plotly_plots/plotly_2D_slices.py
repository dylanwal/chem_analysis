import copy
from typing import Sequence

import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.plotting.plotly_plots.plotly_peaks import plotly_peaks
from chem_analysis.analysis.peak_result import ResultPeaks2D


def plotly_slices(
        signal: Signal2D,
        slices: None | int | Sequence[int] | slice,
        fig: go.Figure | None,
        raw: bool,
        plot_kwargs: dict
) -> go.Figure:
    fig = input_check(fig)

    if raw:
        name = signal.name + "_raw"
        x = signal.x_raw
        y = signal.y_raw
        z = signal.z_raw
    else:
        x = signal.x
        y = signal.y
        z = signal.z
        name = signal.name

    indexes = get_index_from_slice(slices, len(y))
    for i in indexes:
        fig.add_scatter(x=x, y=z[i, :], name=f"{name} (y={i})", **plot_kwargs)
    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.z_label)
    return fig


def get_index_from_slice(slices: None | int | Sequence[int] | slice, len_y: int) -> Sequence[int]:
    if slices is None:
        return range(len_y)
    if isinstance(slices, int):
        return [slices]
    if isinstance(slices, slice):
        return range(
            slices.start if slices.start is not None else 0,
            slices.stop if slices.stop is not None else len_y,
            slices.step if slices.step is not None else 1
        )
    if isinstance(slices, Sequence):
        return slices

    raise TypeError(f"Unsupported type for 'slice'."
                    f"\n\tGiven: {type(slices)}"
                    f"\n\tExpected: None | int | Sequence[int] | slice"
                    )


def plotly_slices_separate(
        signal: Signal2D,
        slices: None | int | Sequence[int] | slice,
        fig: go.Figure | None,
        raw: bool,
        plot_kwargs: dict,
) -> list[go.Figure]:
    fig = input_check(fig)

    if raw:
        name = signal.name + "_raw"
        x = signal.x_raw
        y = signal.y_raw
        z = signal.z_raw
    else:
        x = signal.x
        y = signal.y
        z = signal.z
        name = signal.name

    figs = []
    indexes = get_index_from_slice(slices, len(y))
    for i in indexes:
        fig_ = copy.copy(fig)
        fig_.add_scatter(x=x, y=z[i, :], name=f"{name} (y={i})", **plot_kwargs)
        fig_.layout.xaxis.title = bold_in_html(signal.x_label)
        fig_.layout.yaxis.title = bold_in_html(signal.z_label)
        figs.append(fig_)

    return figs


def plotly_2D_slices_peaks(
        peaks: ResultPeaks2D,
        slices: None | int | Sequence[int] | slice,
        fig: go.Figure | Sequence[go.Figure] | None,
        mode: int | Sequence[int],
        normalize: int,
        plot_kwargs: dict,
) -> go.Figure | list[go.Figure]:
    if fig is None:
        fig = go.Figure()

    indexes = get_index_from_slice(slices, len(peaks.results))
    if isinstance(fig, go.Figure):
        for i in indexes:
            plotly_peaks(peaks.results[i], fig=fig, mode=mode,
                         normalize=normalize, plot_kwargs=copy.copy(plot_kwargs))

        return fig

    if isinstance(fig, Sequence) and isinstance(fig[0], go.Figure) and len(indexes) == len(fig):
        figs = []
        for i, index in enumerate(indexes):
            figs.append(
                plotly_peaks(peaks.results[index], fig=fig[i], mode=mode,
                             normalize=normalize, plot_kwargs=copy.copy(plot_kwargs))
            )
        return figs

    raise TypeError("Invalid input.")
