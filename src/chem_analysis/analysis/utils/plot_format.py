from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go


layout_figure = {
        "autosize": False,
        "width": 1200,
        "height": 600,
        "font": dict(family="Arial", size=18, color="black"),
        "plot_bgcolor": "white"
    }

layout_xaxis = {
    # "title": "<b>X<b>",
    "tickprefix": "<b>",
    "ticksuffix": "</b>",
    "showline": True,
    "linewidth": 5,
    "mirror": True,
    "linecolor": 'black',
    "ticks": "outside",
    "tickwidth": 4,
    "showgrid": False,
    "gridwidth": 1,
    "gridcolor": 'lightgray'
}

layout_yaxis = {
    # "title": "<b>Y<b>",
    "tickprefix": "<b>",
    "ticksuffix": "</b>",
    "showline": True,
    "linewidth": 5,
    "mirror": True,
    "linecolor": 'black',
    "ticks": "outside",
    "tickwidth": 4,
    "showgrid": False,
    "gridwidth": 1,
    "gridcolor": 'lightgray'
}

hex_options = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']


def get_plot_color(num_colors: int) -> Union[list[str], str]:
    """Gets color for 2D plots."""
    color_list = [
        'rgb(10,36,204)',  # blue
        'rgb(172,24,25)',  # red
        'rgb(6,127,16)',  # green
        'rgb(251,118,35)',  # orange
        'rgb(145,0,184)',  # purple
        'rgb(255,192,0)'  # yellow
    ]
    if num_colors <= 1:
        return color_list[0]
    if num_colors <= len(color_list):
        return color_list[:num_colors]
    else:
        num_colors_extra = num_colors - len(color_list)
        for i in range(num_colors_extra):
            color = ["#" + ''.join([np.random.choice(hex_options) for _ in range(6)])]
            color = [col.lstrip('#') for col in color]
            color = ["rgb" + str(tuple(int(col[i:i + 2], 16) for i in (0, 2, 4))) for col in color]
            color_list = color_list + color

        return color_list


#######################################################################################################################
########################### Pandas ####################################################################################
def plot(df: pd.DataFrame, fig=None, auto_open: bool = True, auto_format: bool = True):
    """ Basic plotting """
    if fig is None:
        fig = go.Figure()

    # generate lines
    colors = get_plot_color(len(df.columns))
    traces = get_traces(df, colors)
    for trace in traces:
        fig.add_trace(trace)

    if auto_format:
        add_plot_format(fig, df.index.name, df.columns[0])

    if auto_open:
        fig.write_html("temp.html", auto_open=True)


def get_traces(df: pd.DataFrame, colors: list[str], y_axis_labels: dict = None):
    traces = []
    for i, col in enumerate(df.columns):
        plot_kwargs = {
            "x": df.index,
            "y": df[col],
            "mode": 'lines',
            "connectgaps": True,
            "name": col,
            "line": dict(color=colors[i])
        }
        if y_axis_labels is not None:
            for label in y_axis_labels:
                if col.startswith(label):
                    plot_kwargs["yaxis"] = y_axis_labels[label]
                    break
            else:
                plot_kwargs["yaxis"] = f"y{len(y_axis_labels) + 1}"
                break

        traces.append(go.Scatter(**plot_kwargs))

    return traces


def add_plot_format(fig: go.Figure, x_axis: str, y_axis: str):
    """ Add default format to plot. """
    fig.update_layout(layout_figure)

    # x-axis
    x_axis_format = {"title": x_axis}
    x_axis_format = {**x_axis_format, **layout_xaxis}
    fig.update_xaxes(x_axis_format)

    # y-axis
    y_axis_format = {"title": y_axis}
    y_axis_format = {**y_axis_format, **layout_yaxis}
    fig.update_yaxes(y_axis_format)


def plot_sep_y(df: pd.DataFrame, fig=None, separate_y_axis=None,
               auto_open: bool = True, auto_format: bool = True, spread=0.20):

    """ Basic plotting """
    if fig is None:
        fig = go.Figure()

    # generate lines
    colors = get_plot_color(len(df.columns))
    traces = get_traces(df, colors, separate_y_axis)
    for trace in traces:
        fig.add_trace(trace)

    # adding multiple y-axis
    axis_format = get_multi_y_axis(colors, fig)
    fig.update_layout(**axis_format)

    if auto_open:
        fig.write_html("temp.html", auto_open=True)

    return fig


def get_multi_y_axis(colors: list[str], fig: go.Figure, spread: float = 0.2) -> dict:
    y_axis_labels = {data.yaxis for data in fig.data}
    num_y_axis = len(y_axis_labels)

    gap = spread / (num_y_axis - 1)
    # first trace
    axis_format = {
        "xaxis": {"domain": [spread, 1]},  # 'tickformat': '%I:%M %p \n %b %d'
        # https://github.com/d3/d3-time-format/tree/v2.2.3#locale_format
        "yaxis1": {
            "title": {"text": f"{y_axis_labels[0]}", "standoff": 0},
            # "tickformat": ".4f",
            # "titlefont": {"color": colors[0]},
            # "tickfont": {"color": colors[0]},
            "tickangle": -45
        }
    }
    # the rest of the traces
    for i in range(1, num_y_axis):
        axis_format[f"yaxis{i + 1}"] = {
            "title": {"text": f"{y_axis_labels[i]}", "standoff": 0},
            # "tickformat": ".4f",
            # "titlefont": {"color": colors[i]},
            # "tickfont": {"color": colors[i]},
            "anchor": "free",
            "overlaying": "y",
            "side": "left",
            "position": spread - gap * i,
            "tickangle": -45
        }

    return axis_format


#######################################################################################################################
########################### Pandas ####################################################################################

def plot_series(data: pd.Series, fig=None, auto_open: bool = True, auto_format: bool = True):
    """ Basic plotting """
    if fig is None:
        fig = go.Figure()

    colors = get_plot_color(1)
    fig.add_trace(
        go.Scatter(x=data.index, y=data, mode="lines", connectgaps=True, name=data.name, line=dict(color=colors)))

    if auto_format:
        add_plot_format(fig, data.index.name, data.name)

    if auto_open:
        fig.write_html(f'temp.html', auto_open=True)

    return fig
