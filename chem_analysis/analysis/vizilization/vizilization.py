def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
         op_peaks: bool = True, y_label: str = None, title: str = None, **kwargs) -> go.Figure:
    """ Plot

    General plotting

    Parameters
    ----------
    fig: go.Figure
        plotly figure; will automatically create if not provided
    auto_open: bool
        create "temp.html" and auto_open in browser
    auto_format: bool
        apply built-in formatting
    op_peaks: bool
        add peak plotting stuff
    y_label: str
        y_axis label (used for multiple y-axis)
    title: str
        title

    Returns
    -------
    fig: go.Figure
        plotly figure

    """
    if fig is None:
        fig = go.Figure()

    if "color" in kwargs:
        color = kwargs.pop("color")
    else:
        color = 'rgb(0,0,0)'

    group = self.name

    # add peaks
    if op_peaks:
        if len(self.peaks) > 0:
            if color == 'rgb(0,0,0)':
                peak_color = get_plot_color(self.num_peaks)
            else:
                peak_color = get_similar_color(color, self.num_peaks)

            for peak, color_ in zip(self.peaks, peak_color):
                peak.plot_add_on(fig, color=color_, group=group, y_label=y_label)

    # add main trace
    plot_kwargs = {
        "x": self.result.index,
        "y": self.result,
        "mode": 'lines',
        "connectgaps": True,
        "name": self.result.name,
        "legendgroup": group,
        "line": dict(color=color)
    }
    if y_label is not None:
        plot_kwargs["yaxis"] = y_label

    fig.add_trace(go.Scatter(**plot_kwargs))

    if auto_format:
        if title is not None:
            fig.update_layout(title=title)
        add_plot_format(fig, self.result.index.name, str(self.result.name))

    if auto_open:
        global fig_count
        fig.write_html(f'temp{fig_count}.html', auto_open=True)
        fig_count += 1

    return fig



##################

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             y_label: str = None, title: str = None, **kwargs) -> go.Figure:
        """ Plot

        General plotting

        Parameters
        ----------
        fig: go.Figure
            plotly figure; will automatically create if not provided
        auto_open: bool
            create "temp.html" and auto_open in browser
        auto_format: bool
            apply built-in formatting
        y_label: str
            y_axis label (used for multiple y-axis)
        title: str
            title

        Returns
        -------
        fig: go.Figure
            plotly figure

        """
        if fig is None:
            fig = go.Figure()

        if "color" in kwargs:
            color = kwargs.pop("color")
        else:
            color = 'rgb(0,0,0)'

        # add main trace
        plot_kwargs = {
            "x": x,
            "y": y,
            "mode": 'lines',
            "connectgaps": True,
            "name": self._parent.result.name,
            "line": dict(color=color)
        }
        if y_label is not None:
            plot_kwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(**plot_kwargs))

        if auto_format:
            if title is not None:
                fig.update_layout(title=title)
            add_plot_format(fig, self._parent.result.index.name, str(self._parent.result.name))

        if auto_open:
            global fig_count
            fig.write_html(f'temp{fig_count}.html', auto_open=True)
            fig_count += 1

        return fig

    def plot_add_on(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plot

        Plots several things related to a peak.

        Parameters
        ----------
        fig: go.Figure
            plotly figure
        group: str
            the legend table group that the plot_add_on will be associated with (need to make things appear/disappear together)
        y_label: str
            Label for the y-axis the data is associated with

        """
        self._plot_max(fig, group, y_label, **kwargs)
        # self._plot_bounds(fig, group, **kwargs)
        self._plot_shade(fig, group, y_label, **kwargs)

    def _plot_shade(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plots the shaded area for the peak. """
        kwargs_ = {
            "width": 0
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="lines",
            fill='tozeroy',
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))

    def _plot_max(self, fig: go.Figure, group: str = None, y_label: str = None, **kwargs):
        """ Plots peak name at max. """
        kwargs_ = {
            "size": 1
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(
            x=[self.max_loc],
            y=[self.max],
            mode="text",
            marker=kwargs_,
            text=[f"{self._parent.name}: {self.id_}"],
            textposition="top center",
            showlegend=False,
            **kkwargs
        ))

    def _plot_bounds(self, fig: go.Figure, group: str = None, height: float = 0.08, y_label: str = None, **kwargs):
        """ Adds bounds at the bottom of the plot_add_on for peak area. """
        kwargs_ = {
            "width": 5
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}
        if group:
            kkwargs["legendgroup"] = group
        if y_label:
            kkwargs["yaxis"] = y_label

        bound_height = max(self._parent.result) * height
        # bounds
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.lb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs

        ))
        fig.add_trace(go.Scatter(
            x=[self.hb_loc, self.hb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.hb_loc],
            y=[0, 0],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))
