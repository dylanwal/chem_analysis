
import matplotlib.pyplot as plt


def input_check(fig: plt.Figure | None) -> plt.Figure:
    if fig is None:
        fig = plt.figure()
        return fig
    if not isinstance(fig, plt.Figure):
        raise ValueError("'fig' must be a Matplotlib 'plt.Figure'.")

    return fig
