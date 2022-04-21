import math

import numpy as np
import plotly.graph_objs as go

from chem_analysis.analysis.utils.plot_format import add_plot_format


class PolyLogNorm:
    _count = 0

    def __init__(self, Mn: float = 10_000, D: float = 1.1, Mw: float = None, m_i: np.ndarray = None,
                 name: str = None):
        if 0 >= Mn:
            raise ValueError("Mn must be positive.")
        if Mw is not None and 0 >= Mw:
            raise ValueError("Mw must be positive.")
        if 1 >= D:
            raise ValueError("D must be greater than  1.")

        self.Mn = Mn
        if Mw:
            self.Mw = Mw
            self.D = Mw / Mn
        else:
            self.D = D
            self.Mw = Mn * D

        self.ln_sigma = math.sqrt(math.log(D))
        self.ln_mean = math.log(self.Mn) - self.ln_sigma ** 2 / 2

        if m_i is None:
            m_i = np.logspace(2, 6, 1000)
        self.m_i = m_i

        if name is None:
            self.name = f"log_normal_{self._count}"
            PolyLogNorm._count += 1

    def __repr__(self) -> str:
        return f"Mn: {self.Mn}; D: {self.D}"

    def __call__(self, m_i: np.ndarray = None) -> np.ndarray:
        return self.weight(m_i)

    @property
    def mole(self, m_i: np.ndarray = None) -> np.ndarray:
        if m_i is None:
            m_i = self.m_i
        return 1 / (m_i * np.sqrt(2 * np.pi * np.log(self.D))) * \
               np.exp(-(np.log(m_i / self.Mn) + np.log(self.D) / 2) ** 2 / (2 * np.log(self.D)))

    @property
    def weight(self, m_i: np.ndarray = None) -> np.ndarray:
        if m_i is None:
            m_i = self.m_i
        return 1 / (self.Mn * np.sqrt(2 * np.pi * np.log(self.D))) * \
               np.exp(-(np.log(m_i / self.Mn) + np.log(self.D) / 2) ** 2 / (2 * np.log(self.D)))

    def plot(self, fig: go.Figure = None, op_type: str = "weight", op_log: bool = True,
             auto_open: bool = True, auto_format: bool = True, **kwargs) -> go.Figure:
        """ Basic plotting """
        if fig is None:
            fig = go.Figure()

        if "color" in kwargs:
            color = kwargs.pop("color")
        else:
            color = 'rgb(0,0,0)'

        # add main trace
        plot_kwargs = {
            "mode": 'lines',
            "connectgaps": True,
            "name": self.name,
            "line": dict(color=color)
        }
        if op_type == "weight":
            plot_kwargs["x"] = self.m_i
            plot_kwargs["y"] = self.weight
            y_axis = "weight fraction"
        elif op_type == "mole":
            plot_kwargs["x"] = self.m_i
            plot_kwargs["y"] = self.mole
            y_axis = "mole fraction"
        else:
            raise ValueError(f"Invalid 'op_type'. (given: {op_type})")

        fig.add_trace(go.Scatter(**plot_kwargs))

        if auto_format:
            if op_log:
                format_kwargs = {"type": "log"}
            else:
                format_kwargs = {}
            add_plot_format(fig, "molecular weight", y_axis, x_kwargs=format_kwargs)

        if auto_open:
            fig.write_html(f'temp.html', auto_open=True)

        return fig


def run_local():
    dis = PolyLogNorm()
    dis.plot()


if __name__ == '__main__':
    run_local()
