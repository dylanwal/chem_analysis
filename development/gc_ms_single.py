import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)
INTERNAL_STANDARD = LIBRARY.find_by_label("TCB")
picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")


def main(data_path):
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(data_path)

    # fid
    fid.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.7), invert=True))
    fid.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(fid,
                                                               scipy_kwargs={"height": 8000, "width": 0.1}
                                                               )
    fid_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                     n_points_with_pos_slope=2)
    fid_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_fid, fid_peaks)

    fid_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(fid, fig=fid_fig)
    ca.plotting.peaks(fid_compounds, fig=fid_fig)
    fid_fig.layout.title = "FID"

    # ms
    ms.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(ms,
                                                               scipy_kwargs={"prominence": 10000}
                                                               )
    ms_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                    n_points_with_pos_slope=2)
    ms_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_ms, ms_peaks)

    # plotting peak results
    ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=ms_fig)
    ca.plotting.peaks(ms_compounds, fig=ms_fig)
    ms_fig.layout.title = "MS"

    ca.plotting.PlotlyConfig.merge_figures([fid_fig, ms_fig])


if __name__ == "__main__":
    data_path_ = r"C:\Users\nicep\Desktop\11_23\DJW-11-23-90min-TMS.D"
    main(data_path_)
