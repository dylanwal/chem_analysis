
import numpy as np
import plotly.graph_objs as go

from chem_analysis.plotting.plotly_plots.plotly_utils import input_check
from chem_analysis.sec.sec_calibration import SECCalibration


def plotly_calibration(calibration: SECCalibration, plot_kwargs: dict, fig: go.Figure | None):
    fig = input_check(fig)
    time_ = np.linspace(calibration.x_bounds[0], calibration.x_bounds[1], 100)
    mw = calibration.get_y(time_)

    kwargs_ = dict(
        name="cal",
        mode="lines",
        yaxis="y2",
        legendgroup="cal",
    )
    plot_kwargs = kwargs_ | plot_kwargs
    plot_kwargs["line"] = {"width": 1, "color": "rgb(130,130,130)"} | plot_kwargs.get("line", dict())

    fig.add_trace(go.Scatter(
        x=time_,
        y=mw,
        showlegend=True,
        **plot_kwargs,
    ))

    plot_kwargs["line"] = {"dash": 'dash'} | plot_kwargs["line"]

    # low limit
    fig.add_trace(go.Scatter(
            x=[calibration.x_bounds[0], calibration.x_bounds[0]],
            y=[0, np.max(mw)],
            showlegend=False,
            **plot_kwargs,
        ))

    # up limit
    fig.add_trace(go.Scatter(
            x=[calibration.x_bounds[1], calibration.x_bounds[1]],
            y=[0, np.max(mw)],
            showlegend=False,
            **plot_kwargs,
        ))

    fig.update_layout(
        yaxis2=dict(
            title="molecular weight",
            titlefont=dict(
                color=plot_kwargs["line"]["color"],
            ),
            tickfont=dict(
                color=plot_kwargs["line"]["color"],
            ),
            anchor="x",
            overlaying="y",
            side="right",
            type="log",
            range=[2, 6]
        ),
    )
    return fig
