
import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_2d import Signal2D


def plotly_surface(
        signal: Signal2D,
        fig: go.Figure | None,
        raw: bool,
        plot_kwargs: dict,
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

    fig.add_surface(x=x, y=y, z=z, name=name, **plot_kwargs)
    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    fig.layout.yaxis.title = bold_in_html(signal.z_label)
    return fig
