import pathlib
from typing import Iterable

import plotly.graph_objs as go
import plotly.colors as colors
import plotly.io as pio


def plotly_template():
    """ use: fig = go.Figure(layout=dict(template=plotly_template())) """
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


pio.templates.default = plotly_template()


def plotly_layout() -> dict:
    return {"template": plotly_template()}


def color_scales_figure() -> go.Figure:
    import plotly.express as px
    return px.colors.sequential.swatches_continuous()


def get_colors_from_scale(n: int = 10, color_scale: str = 'Viridis') -> list[str]:
    color_scale = colors.get_colorscale(color_scale)
    return colors.sample_colorscale(color_scale, [i / (n - 1) for i in range(n)])


def merge_html_figs(
        figs: list[go.Figure | str],
        filename: str | pathlib.Path = "merged_htmls.html",
        auto_open: bool = True,
        title: str | None = None,
        html_head: str | Iterable[str] | None = None,
):
    """
    Merges plotly figures into single html

    Parameters
    ----------
    figs: list[go.Figure, str]
        list of figures to append together or html divs
    filename: str
        file name
    auto_open: bool
        open html in browser after creating
    title: str | None
        title of the figure
    html_head: str | Iterable[str] | None
        headers to add to html
    """
    if not isinstance(filename, pathlib.Path):
        filename = pathlib.Path(filename)
    if filename.suffix != ".html":
        filename = filename.with_suffix(".html")

    head = '\n\t<meta charset="UTF-8">\n\t<meta name="viewport" content="width=device-width, initial-scale=1.0">'
    if title is not None:
        head += f"\n\t<title>{title}</title>"
    if html_head is not None:
        if isinstance(html_head, str):
            html_head = [html_head]
        for header in html_head:
            head += f"\n\t{header}"

    body = ""
    if title is not None:
        body += f"\n\t<h1>{title}</h1>"
    for fig in figs:
        if isinstance(fig, str):
            body += "\n\t" + fig
            continue

        inner_html = fig.to_html(include_plotlyjs="cdn").split('<body>')[1].split('</body>')[0]
        body += inner_html

    text = f'<!DOCTYPE html>\n<html lang="en">\n<head>{head}</head>\n<body>{body}</body>'

    with open(filename, 'w') as file:
        file.write(text)

    if auto_open:
        import os
        os.system(fr"start {filename}")