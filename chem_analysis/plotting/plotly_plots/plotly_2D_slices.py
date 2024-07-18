from typing import Sequence

import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_2d import Signal2D


def plotly_slices(
        signal: Signal2D,
        slices: None | int | Sequence[int] | slice,
        fig: go.Figure | None,
        config: PlotlyConfig | None = None,
        raw: bool = True,
) -> go.Figure:
    fig, config = PlotlyConfig.input_check(fig, config)

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
        fig.add_scatter(x=x, y=z[i, :], name=f"{name} (y={i})")
    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    fig.layout.yaxis.title = bold_in_html(signal.z_label)
    return fig


def get_index_from_slice(slices: None | int | Sequence[int] | slice, len_y: int) -> Sequence[int]:
    if slices is None:
        return range(len_y)
    if isinstance(slices, int):
        return [slices]
    if isinstance(slices, slice):
        return range(slices.start, slices.stop)
    if isinstance(slices, Sequence):
        return slices

    raise TypeError(f"Unsupported type for 'slice'."
                    f"\n\tGiven: {type(slices)}"
                    f"\n\tExpected: None | int | Sequence[int] | slice"
                    )
