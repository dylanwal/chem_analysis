
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
from chem_analysis.processing.processing_method import Baseline


def plotly_baseline(baseline: Baseline, fig: go.Figure | None, config: PlotlyConfig | None) -> go.Figure:
    fig, config = PlotlyConfig.input_check(fig, config)

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.data,
            mode="lines",
            name="raw_signal",
            connectgaps=config.signal_connect_gaps,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.baseline,
            mode="lines",
            name="baseline",
            connectgaps=config.signal_connect_gaps,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=baseline.x,
            y=baseline.data - baseline.baseline,
            mode="lines",
            name="result",
            connectgaps=config.signal_connect_gaps,
        )
    )

    return fig


