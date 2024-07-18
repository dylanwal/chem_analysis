
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.processing.processing_method import Baseline


def plotly_baseline(
        baseline: Baseline,
        plot_kwargs: dict,
        fig: go.Figure | None,
) -> go.Figure:
    fig = input_check(fig)

    kwargs = dict(
            mode="lines",
            connectgaps=True,
    )
    plot_kwargs = kwargs | plot_kwargs

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.data,
            name="raw_signal",
            **plot_kwargs
        )
    )

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.baseline,
            name="baseline",
            **plot_kwargs
        )
    )

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.data - baseline.baseline,
            name="result",
            **plot_kwargs
        )
    )

    return fig


