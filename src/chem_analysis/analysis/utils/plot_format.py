from typing import Union

import numpy as np


layout_figure = {
        "autosize": False,
        "width": 900,
        "height": 600,
        # "showlegend": False,
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


def get_plot_color(num_colors: int) -> Union[list[str], str]:
    """Gets color for 2D plots."""
    color_list = [
        'rgb(10, 36, 204)',  # blue
        'rgb(172, 24, 25)',  # red
        'rgb(6, 127, 16)',  # green
        'rgb(251, 118, 35)',  # orange
        'rgb(145, 0, 184)',  # purple
        'rgb(255, 192, 0)'  # yellow
    ]
    if num_colors <= 1:
        return color_list[0]
    if num_colors <= len(color_list):
        return color_list[:num_colors]
    else:
        num_colors_extra = num_colors - len(color_list)
        for i in range(num_colors_extra):
            color = ["#" + ''.join([np.random.choice('0123456789ABCDEF') for _ in range(6)])]
            color = [col.lstrip('#') for col in color]
            color = ["rgb" + str(tuple(int(col[i:i + 2], 16) for i in (0, 2, 4))) for col in color]
            color_list = color_list + color

        return color_list
