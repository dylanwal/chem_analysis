from __future__ import annotations

import plotly.graph_objs as go

from chem_analysis.plotting.config import PlotConfig


class PlotlyConfig(PlotConfig):
    def __init__(self):
        super().__init__()

    @classmethod
    def input_check(cls, fig: go.Figure | None, config: PlotlyConfig | None) -> tuple[go.Figure, PlotlyConfig]:
        if fig is None:
            fig = go.Figure()
        else:
            if not isinstance(fig, go.Figure):
                raise ValueError("'fig' must be a plotly 'go.Figure'.")
        if config is None:
            config = PlotConfig()
        else:
            if not isinstance(config, PlotlyConfig):
                raise ValueError("'config' must be a 'PlotlyConfig'.")

        return fig, config

    @staticmethod
    def plotly_layout():
        import plotly.graph_objs as go
        template = go.layout.Template()
        template.layout.font = dict(family="Arial", size=18, color="black")
        template.layout.plot_bgcolor = "white"
        template.layout.width, template.layout.height = 1200, 600
        template.layout.xaxis.tickprefix = "<b>"
        template.layout.xaxis.ticksuffix = "<b>"
        template.layout.xaxis.showline = True
        template.layout.xaxis.linewidth = 5
        template.layout.xaxis.linecolor = "black"
        template.layout.xaxis.ticks = "outside"
        template.layout.xaxis.tickwidth = 4
        template.layout.xaxis.showgrid = False
        template.layout.xaxis.mirror = True
        template.layout.yaxis.tickprefix = "<b>"
        template.layout.yaxis.ticksuffix = "<b>"
        template.layout.yaxis.showline = True
        template.layout.yaxis.linewidth = 5
        template.layout.yaxis.linecolor = "black"
        template.layout.yaxis.ticks = "outside"
        template.layout.yaxis.tickwidth = 4
        template.layout.yaxis.showgrid = False
        template.layout.yaxis.mirror = True

        return template
