import time
import json

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


def generate_data(n: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(0, 100, n)
    signal = np.exp(-0.1 * (x - 80) ** 2) + np.exp(-2 * (x - 20) ** 2) + 0.05 * x + 0.2 * np.sin(0.1 * x) - 3*np.exp(-0.1 * (x - 50) ** 2)  # Peak + Drift
    noise = 0.05 * np.random.randn(len(x))  # Add noise
    return x, signal + noise


def generate_data2(n: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(0, 100, n)
    signal = np.exp(-0.1 * (x - 80) ** 2) + np.exp(-2 * (x - 20) ** 2) + 0.05 * x + 0.2 * np.sin(0.1 * x) + 3*np.exp(-0.1 * (x - 50) ** 2)  # Peak + Drift
    noise = 0.05 * np.random.randn(len(x))  # Add noise
    return x, signal + noise


def create_figure(x: np.ndarray, y: np.ndarray) -> tuple[go.Figure, dict]:
    fig = go.Figure(layout={'template': 'simple_white'})
    fig.add_scatter(x=x, y=y, name='raw_data', mode='lines', line={'width': 5})

    props = dict()

    start = time.perf_counter()
    method = ca.p.baseline.Polynomial(degree=3, save_result=True)
    method.run_xy(x, y)
    end = time.perf_counter()
    props['Polynomial'] = {'time': end - start}
    fig.add_scatter(x=x, y=method.baseline, name='Polynomial', mode='lines')

    start = time.perf_counter()
    method = ca.p.baseline.MorphologicalAverage(save_result=True)
    method.run_xy(x, y)
    end = time.perf_counter()
    props['MorphologicalAverage'] = {'time': end - start}
    fig.add_scatter(x=x, y=method.baseline, name='MorphologicalAverage', mode='lines')

    start = time.perf_counter()
    method = ca.p.baseline.SectionMinMax(save_result=True)
    method.run_xy(x, y)
    end = time.perf_counter()
    props['SectionMinMax'] = {'time': end - start}
    fig.add_scatter(x=x, y=method.baseline, name='SectionMinMax', mode='lines')

    start = time.perf_counter()
    method = ca.p.baseline.RLOESS(save_result=True)
    method.run_xy(x, y)
    end = time.perf_counter()
    props['RLOESS'] = {'time': end - start}
    fig.add_scatter(x=x, y=method.baseline, name='RLOESS', mode='lines')

    start = time.perf_counter()
    method = ca.p.baseline.Wavelet(save_result=True)
    method.run_xy(x, y)
    end = time.perf_counter()
    props['Wavelet'] = {'time': end - start}
    fig.add_scatter(x=x, y=method.baseline, name='Wavelet', mode='lines')

    return fig, props


def main():
    n = 1000
    figs = []

    x, y = generate_data(n)
    fig, props = create_figure(x, y)
    figs.append(fig)
    figs.append("<p>" + json.dumps(props, indent=2) + "</p>")

    x, y = generate_data2(n)
    fig, props = create_figure(x, y)
    figs.append(fig)
    figs.append("<p>" + json.dumps(props, indent=2) + "</p>")

    import chem_analysis.plotting.plotly_plots.plotly_utils as plotly_utils
    plotly_utils.merge_figures(figs, filename=r"./figs/baseline.html", auto_open=True)


if __name__ == '__main__':
    main()
