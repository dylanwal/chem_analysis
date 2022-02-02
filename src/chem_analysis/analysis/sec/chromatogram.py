from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from Dylan.sec.plot_format import layout_figure, layout_xaxis, layout_yaxis
from Dylan.sec.trace import Trace


class Chromatogram:
    """
    traces need
    """
    def __init__(self, df: pd.DataFrame):
        self.df = df

    def plot(self, fig=None, auto_open: bool = True, auto_format: bool = True, sep_y_axis: Union[bool, list] = False,
              kwargs_lines: dict = None, kwargs_layout: dict = None, kwargs_x_axes: dict = None, kwargs_y_axes: dict =
             None):
        """ Basic plotting """
        if fig is None:
            fig = go.Figure()

        # generate lines
        traces = self._get_traces(kwargs_lines)
        for trace in traces:
            fig.add_trace(trace)

        #  add format
        if auto_format:
            self._plot_format(fig)



            fig.add_trace(go.Scatter(**plot_kwargs))

        if separate_y_axis:
            spread = .20
            gap = spread / (len(columns) - 1)
            axis_format = {
                "xaxis": {"domain": [spread, 1], 'tickformat': '%I:%M %p \n %b %d'},
                # https://github.com/d3/d3-time-format/tree/v2.2.3#locale_format
                "yaxis1": {
                    "title": {"text": f"{columns[0]}", "standoff": 0},
                    "titlefont": {"color": colors[0]},
                    "tickfont": {"color": colors[0]},
                    "tickangle": -45
                }
            }
            for i in range(1, len(columns)):
                axis_format[f"yaxis{i + 1}"] = {
                    "title": {"text": f"{columns[i]}", "standoff": 0},
                    "titlefont": {"color": colors[i]},
                    "tickfont": {"color": colors[i]},
                    "anchor": "free",
                    "overlaying": "y",
                    "side": "left",
                    "position": spread - gap * i,
                    "tickangle": -45
                }

            fig.update_layout(**axis_format)

        if auto_open:
            fig.write_html("temp.html", auto_open=True)

        return fig

    def _get_traces(self):
        if kwargs:
            if "line" in kwargs:
                if "color" not in kwargs["line"]:
                    kwargs["line"] = {"color": get_plot_color(1)}
            else:
                kwargs["line"] = {"color": get_plot_color(1)}
        else:
            kwargs = {"line": {"color": get_plot_color(1)}}

        colors = get_color_for_plotting(len(columns))
        for i, col in enumerate(columns):
            plot_kwargs = {
                "x": data[x_axis],
                "y": data[col],
                "mode": 'lines',
                "connectgaps": True,
                "name": col,
                "line": dict(color=colors[i])
            }
            if separate_y_axis:
                plot_kwargs["yaxis"] = f"y{i + 1}"


    def _plot_format(self, fig):
        """ Add default format to plot. """
        fig.update_layout(layout_figure)

        # x-axis
        if self.x_unit is not None:
            x_axis_format = {"title": f"<b>{self.x_label}<b> ({self.x_unit})"}
        else:
            x_axis_format = {"title": self.x_label}
        x_axis_format = {**x_axis_format, **layout_xaxis}
        fig.update_xaxes(x_axis_format)

        # y-axis
        if self.y_unit is not None:
            y_axis_format = {"title": f"<b>{self.y_label}<b> ({self.y_unit})"}
        else:
            y_axis_format = {"title": self.y_label}
        y_axis_format = {**y_axis_format, **layout_yaxis}
        fig.update_yaxes(y_axis_format)

    def plot_2D_df(data: pd.DataFrame, x_axis: str = None, fig: go.Figure = None, layout_kwargs: dict = None,
                   xaxes_kwargs: dict = None, yaxes_kwargs: dict = None, separate_y_axis: bool = False,
                   auto_open: bool =
                   True) -> go.Figure:

        if fig is None:
            fig = go.Figure()


