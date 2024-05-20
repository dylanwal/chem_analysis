
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)


def main():
    chem = LIBRARY.find_by_label("TCB")
    ms = chem.get_ms()

    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=fig)
    fig.layout.title = f"{chem.name}"
    fig.write_html(f"ms_{chem.label}.html", include_plotlyjs='cdn')


if __name__ == "__main__":
    main()
