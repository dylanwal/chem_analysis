import pathlib

import plotly.graph_objs as go

import chem_analysis as ca


def main():
    file_path = pathlib.Path(r"data//ms.csv")
    ms = ca.ms.MSSignal.from_csv(file_path)

    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=fig)
    fig.layout.title = f"ms of 1,2,3-Trichlorobenzene"
    fig.show()
    fig.write_image('figs/ms_data_analysis.png', width=fig.layout.template.layout.width, height=fig.layout.template.layout.height)


if __name__ == "__main__":
    main()
