
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.gc_lc.GCLibrary.from_JSON(lib_path)
# LIBRARY = ca.gc_lc.GCLibrary.from_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\nist20_1.json")


def main():
    chem = LIBRARY.find_by_label("TCB")
    print(chem.methods)
    # ms = chem.get_ms()
    #
    # fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # ca.plotting.signal(ms, fig=fig)
    # fig.layout.title = f"{chem.name}"
    # fig.write_html(f"ms_{chem.label}.html", include_plotlyjs='cdn')
    print("hi")


if __name__ == "__main__":
    main()
