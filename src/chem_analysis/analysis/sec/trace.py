
import numpy as np
import plotly.graph_objs as go

from Dylan.sec.src.sec.utils.plot_format import layout_figure, layout_xaxis, layout_yaxis, get_plot_color


class Trace:
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 x_unit: str = None,
                 y_unit: str = None,
                 **kwargs
                 ):
        self.x = x
        self.y = y
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "y_axis"
        self.x_unit = x_unit
        self.y_unit = y_unit

        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        text = ""
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self.x)})"
        return text

    def plot(self, fig=None, auto_open: bool = True, auto_format: bool = True, **kwargs):
        """ Basic plotting """
        if fig is None:
            fig = go.Figure()

        if kwargs:
            if "line" in kwargs:
                if "color" not in kwargs["line"]:
                    kwargs["line"] = {"color": get_plot_color(1)}
            else:
                kwargs["line"] = {"color": get_plot_color(1)}
        else:
            kwargs = {"line": {"color": get_plot_color(1)}}

        fig.add_trace(
            go.Scatter(x=self.x, y=self.y, mode="lines", **kwargs)
        )

        if auto_format:
            self._plot_format(fig)

        if auto_open:
            fig.write_html(f'temp.html', auto_open=True)

        return fig

    def _plot_format(self, fig):
        """ Add default format to plot. """
        fig.update_layout(layout_figure)

        # x-axis
        if self.x_unit is not None:
            x_axis_format = {"title": f"<b>{self.x_label}<b> ({self.x_unit})"}
        else:
            x_axis_format = {"title": self.x_label}
        x_axis_format = {**x_axis_format, **layout_xaxis}
        fig.update_xaxes(x_axis_format)

        # y-axis
        if self.y_unit is not None:
            y_axis_format = {"title": f"<b>{self.y_label}<b> ({self.y_unit})"}
        else:
            y_axis_format = {"title": self.y_label}
        y_axis_format = {**y_axis_format, **layout_yaxis}
        fig.update_yaxes(y_axis_format)


def local_run():
    trace = Trace(
        x=np.linspace(0, 100, 100),
        y=np.linspace(0, 100, 100) + 30 * np.random.random(100),
        x_label="x_test",
        y_label="y_test",
        x_unit="per",
        y_unit="mm"
    )
    trace.plot()


if __name__ == '__main__':
    local_run()
