from typing import Any

import chem_analysis.plotting.plotting_primitives as primitives

try:
    import plotly.graph_objs as go
except ImportError:
    raise ImportError(
        "Please install matplotlib with `pip install plotly` and `pip install kaleido==0.1.0post1` "
        "(for image generation) or use another plotting package."
    )

PEAK_COUNTER = 0

def draw_container(fig: go.Figure, container: primitives.DrawingContainer, label_mode: int = 0) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    kwargs = dict(hoverinfo="skip")
    if label_mode == 1:
        global PEAK_COUNTER
        kwargs['showlegend'] = True
        kwargs['legendgroup'] = f"peaks_{PEAK_COUNTER}"
        PEAK_COUNTER += 1
    else:
        kwargs['showlegend'] = False

    draw_containers(fig, container, kwargs)
    return fig


def draw_containers(fig: go.Figure, container: primitives.DrawingContainer, kwargs: dict[str, Any] = None):
    draw_dots(fig, container.dots, kwargs)
    draw_lines(fig, container.lines, kwargs)
    draw_fills(fig, container.fills, kwargs)
    # draw_texts(fig, container.texts, kwargs)


def draw_dots(fig: go.Figure, dots: list[primitives.Dots], kwargs: dict[str, Any] = None):
    for d in dots:
        fig.add_scatter(
            x=d.x,
            y=d.y,
            mode="markers",
            marker=dict(color=d.color, size=d.size),
            **(kwargs or {})
        )


def draw_lines(fig: go.Figure, lines: list[primitives.Lines], kwargs: dict[str, Any] = None):
    for l in lines:
        fig.add_scatter(
            x=l.x,
            y=l.y,
            mode="lines",
            line=dict(color=l.color, width=l.width, dash=l.dash),
            **(kwargs or {})
        )


def draw_fills(fig: go.Figure, fills: list[primitives.Fills], kwargs: dict[str, Any] = None):
    for f in fills:
        fig.add_scatter(
            x=f.x,
            y=f.y,
            mode="lines",
            fill="toself",
            fillcolor=f.color,
            line=dict(width=0),
            **(kwargs or {})
        )


# def draw_texts(fig: go.Figure, text: list[primitives.Texts]):
#     for t in text:
#         t.prepare_for_drawing("\n")
#         fig.add_scatter(
#             x=t.x,
#             y=t.y,
#             text=t.symbols,
#             mode="text",
#             textfont=dict(color=t.color, family=t.font, size=t.size),
#         )