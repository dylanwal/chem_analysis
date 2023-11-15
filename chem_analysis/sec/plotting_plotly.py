
import numpy as np
import plotly.graph_objs as go

from chem_analysis.config import global_config
from chem_analysis.sec.sec_calibration import SECCalibration


def plotly_sec_calibration(calibration: SECCalibration, *, fig: go.Figure, config):
    if fig is None:
        fig = go.Figure()

    # fig.add_trace(go.Scatter(
    #     x=[calibration.x_bounds[0], calibration.x_bounds[0]],
    #     y=[0, 1]
    # ))
    # fig.add_trace(go.Scatter(
    #     x=[calibration.x_bounds[1], calibration.x_bounds[1]],
    #     y=[0, 1]
    # ))
    plotly_sec_calibration_secondary(fig, calibration)

    return fig


def plotly_sec_calibration_secondary(fig: go.Figure, calibration: SECCalibration):
    color = "rgb(130,130,130)"
    kwargs_ = {
        "width": 1,
        "color": color
    }
    time_ = np.linspace(calibration.x_bounds[0], calibration.x_bounds[1], 100)
    mw = calibration.get_y(time_)
    fig.add_trace(go.Scatter(
        x=time_,
        y=mw,
        name="cal",
        mode="lines",
        line=kwargs_,
        yaxis="y2",
        showlegend=True,
        legendgroup="cal",
    ))

    kwargs_["dash"] = 'dash'

    # low limit
    fig.add_trace(go.Scatter(
            x=[calibration.x_bounds[0], calibration.x_bounds[0]],
            y=[0, np.max(mw)],
            name="calibration",
            mode="lines",
            line=kwargs_,
            yaxis="y2",
            showlegend=False,
            legendgroup="cal",
        ))

    # up limit
    fig.add_trace(go.Scatter(
            x=[calibration.x_bounds[1], calibration.x_bounds[1]],
            y=[0, np.max(mw)],
            name="calibration",
            mode="lines",
            line=kwargs_,
            yaxis="y2",
            showlegend=False,
            legendgroup="cal",
        ))

    fig.update_layout(
        yaxis2=dict(
            title="molecular weight",
            titlefont=dict(
                color=color
            ),
            tickfont=dict(
                color=color
            ),
            anchor="x",
            overlaying="y",
            side="right",
            type="log",
            range=[2, 6]
        ),
    )
