
import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.plotting.plot_format import bold_in_html
from chem_analysis.base_obj.signal_2d import Signal2D


def plotly_surface(
        signal: Signal2D,
        fig: go.Figure | None,
        config: PlotlyConfig | None = None,
        raw: bool = True
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

    fig.add_surface(x=x, y=y, z=z, name=name)
    fig.layout.xaxis.title = bold_in_html(signal.x_label)
    fig.layout.yaxis.title = bold_in_html(signal.y_label)
    fig.layout.yaxis.title = bold_in_html(signal.z_label)
    return fig
