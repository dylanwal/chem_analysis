def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
         op_peaks: bool = True, op_cal: bool = True, y_label: str = None, title: str = None, **kwargs) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    fig = super().plot(fig, auto_open=False, op_peaks=op_peaks, y_label=y_label)

    if op_cal and self.cal is not None:
        self._plot_cal(fig)

    if auto_format:
        if title is not None:
            fig.update_layout(title=title)
        add_plot_format(fig, self.result.index.name, str(self.result.name))

    if auto_open:
        global fig_count
        fig.write_html(f'temp{fig_count}.html', auto_open=True)
        fig_count += 1

    return fig


def _plot_cal(self, fig: go.Figure, **kwargs):
    color = "rgb(130,130,130)"
    kwargs_ = {
        "width": 1,
        "color": color
    }
    if kwargs:
        kwargs_ = {**kwargs_, **kwargs}

    kkwargs = {}

    mw = self.cal(self.result.index.to_numpy())
    time = self.result.index
    fig.add_trace(go.Scatter(
        x=time,
        y=mw,
        name=self.cal.name if self.cal.name is not None else "calibration",
        mode="lines",
        line=kwargs_,
        yaxis="y2",
        showlegend=True,
        legendgroup="cal",
        **kkwargs
    ))

    kwargs_["dash"] = 'dash'

    # low limit
    if self.cal.lb:
        fig.add_trace(go.Scatter(
            x=[self.cal.lb_loc, self.cal.lb_loc],
            y=[0, np.max(mw)],
            name="calibration",
            mode="lines",
            line=kwargs_,
            yaxis="y2",
            showlegend=False,
            legendgroup="cal",
            **kkwargs
        ))

    # up limit
    if self.cal.ub:
        fig.add_trace(go.Scatter(
            x=[self.cal.ub_loc, self.cal.ub_loc],
            y=[0, np.max(mw)],
            name="calibration",
            mode="lines",
            line=kwargs_,
            yaxis="y2",
            showlegend=False,
            legendgroup="cal",
            **kkwargs
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

def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
         op_peaks: bool = True, op_cal: bool = True, **kwargs) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    colors = get_plot_color(self.num_signals)

    for i, (sig, color) in enumerate(zip(self, colors)):
        kwargs_ = {"color": color}
        if op_cal and i == 0:
            kwargs_["op_cal"] = True
        else:
            kwargs_["op_cal"] = False
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}
        fig = sig.plot_add_on(fig, auto_open=False, auto_format=False, op_peaks=op_peaks, **kwargs_)

    if auto_format:
        add_plot_format(fig, self.x_label, "; ".join(self.y_labels))

    if auto_open:
        global fig_count
        fig.write_html(f'temp{fig_count}.html', auto_open=True)
        fig_count += 1

    return fig

def plot_mw(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True, y_label: str = None,
            title: str = None, spread: float = None, **kwargs) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    colors = get_plot_color(2)

    # add main trace
    plot_kwargs_xi = {
        "x": self.mw_i,
        "y": self.xi,
        "mode": 'lines',
        "connectgaps": True,
        "name": "xi",
        "line": dict(color=colors[0]),
        "yaxis": "y1"
    }

    fig.add_trace(go.Scatter(**plot_kwargs_xi))
    plot_kwargs_wi = {
        "x": self.mw_i,
        "y": self.wi,
        "mode": 'lines',
        "connectgaps": True,
        "name": "wi",
        "line": dict(color=colors[1]),
        "yaxis": "y2"
    }

    fig.add_trace(go.Scatter(**plot_kwargs_wi))

    # adding multiple y-axis
    y_axis_labels = ["y1", "y2"]
    if spread is None:
        spread = 0.05 * len(set(y_axis_labels))
    axis_format = get_multi_y_axis(colors, fig, spread)
    fig.update_layout(**axis_format)

    if auto_open:
        global fig_count
        fig.write_html(f'temp{fig_count}.html', auto_open=True)
        fig_count += 1

    return fig