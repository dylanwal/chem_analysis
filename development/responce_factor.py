import glob
import pathlib

import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)

REFERENCE_DATA = [
     {
            "label": "K6_2",
            "group": "ketone",
            "name": "2-hexanone",
            "cas": "591-78-6",
            "smiles": "CCCCC(=O)C",
            "retention_time": 4.704,
            "response": 4.67,
            "mass_spectrum": None
        },
        {
            "label": "K7_2",
            "group": "ketone",
            "name": "2-heptanone",
            "cas": "110-43-0",
            "smiles": "CCCCCC(=O)C",
            "retention_time": 8.149,
            "response": 3.57,
            "mass_spectrum": None
        },
        {
            "label": "K8_2",
            "group": "ketone",
            "name": "2-octanone",
            "cas": "111-13-7",
            "smiles": "CCCCCCC(=O)C",
            "retention_time": 13.3,
            "response": 2.84,
            "mass_spectrum": None
        },
        {
            "label": "K9_2",
            "group": "ketone",
            "name": "2-nonanone",
            "cas": "821-55-6",
            "smiles": "CCCCCCCC(=O)C",
            "retention_time": 19.517,
            "response": 2.12,
            "mass_spectrum": None
        },

]


def process_one(signal: ca.base_obj.Signal, type_: str) -> ca.analysis.peak_picking.ResultPeaks:
    if type_ == "fid":
        signal.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.3), invert=True))
    signal.processor.add(ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True))

    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(signal, scipy_kwargs={"height": 20000, "width": 0.1})
    peak_result_fid = ca.analysis.boundary_detection.rolling_ball(peak_locations, n=10, min_height=0.002, n_points_with_pos_slope=1)

    fig = go.Figure(layout=dict(template=ca.plotting.PlotlyConfig.plotly_layout()))
    # ca.plotting.baseline(signal, fig=fig)
    ca.plotting.signal(signal, fig=fig)
    ca.plotting.peaks(peak_result_fid, fig=fig)

    current_file_path = pathlib.Path(__file__)
    fig.write_html(current_file_path.parent / f"figs/{type_}_{signal.name}.html")

    return peak_result_fid


def main(root_folder: str, pattern: str):
    specific_folders = glob.glob(root_folder + pattern)
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]

    fid_peak_results = []
    ms_peak_results = []
    for ms, fid in data:
        fid_peak_results.append(process_one(fid, "fid"))
        ms_peak_results.append(process_one(ms, "ms"))

    print(fid_peak_results)
    #
    # data = ca.gc_lc.GCMSSignal2D.from_file(folder_ + r"\data.npz")
    # results = ResultGrouper(time_=data.y)
    # for sig in range(len(data)):
    #     results.add_result(sig, process_single(data.get_signal(sig)))
    # TCB = chemistry_lib.find_by_label("TCB")
    # results.set_calibrate(TCB)


def main_first():
    lib_ = ca.gc_lc.GCLibrary('DJW-oxidation')
    lib_.to_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON")


if __name__ == "__main__":
    # main_first()
    root_folder_ = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern_ = r"\DJW-cal-Kn_2-*.D"
    main(root_folder_, pattern_)
