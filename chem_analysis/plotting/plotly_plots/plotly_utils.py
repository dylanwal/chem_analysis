from __future__ import annotations

import pathlib
from typing import Iterable, Sequence

import plotly.graph_objs as go


def input_check(fig: go.Figure | None) -> go.Figure:
    if fig is None:
        fig = go.Figure()
    else:
        if not (isinstance(fig, go.Figure) or (isinstance(fig, Sequence) and isinstance(fig[0], go.Figure))):
            raise ValueError("'fig' must be a plotly 'go.Figure'.")

    return fig


def layout() -> dict:
    return dict(template=template())


def template() -> go.Template:
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
    # template.layout.hoverlabel.bgcolor = "white"
    # template.layout.hoverlabel.font.size = 12
    template.layout.hoverlabel.font.family = "Arial"

    return template


def merge_figures(
        figs: list[go.Figure | str],
        title: str | None = None,
        html_head: str | Iterable[str] | None = None,
        auto_open: bool = False,
        filename: str | pathlib.Path | None = None,
) -> str:
    """
        Merges plotly figures into single html

        Parameters
        ----------
        figs: list[go.Figure, str]
            list of figures to append together or html divs
        title: str | None
            title of the figure
        html_head: str | Iterable[str] | None
            headers to add to html
        auto_open: bool
            whether to automatically open the html
        filename: str | pathlib.Path | None
            If provided, the html will be saved to this filename
    """
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

        # inner_html = fig.to_html(include_plotlyjs="cdn").split('<body>')[1].split('</body>')[0]
        inner_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
        body += inner_html

    html = f'<!DOCTYPE html>\n<html lang="en">\n<head>{head}</head>\n<body>{body}</body>'

    if auto_open and filename is None:
        filename = "merged_figs.html"
    if filename is not None:
        if not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)
        if filename.suffix != ".html":
            filename = filename.with_suffix(".html")

        with open(filename, 'w') as file:
            file.write(html)

    if auto_open:
        import os
        os.system(fr"start {filename}")

    return html
