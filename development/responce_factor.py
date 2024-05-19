import glob
import pathlib
import shutil

import numpy as np
from scipy import stats
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)


def process_one_signal(signal: ca.base_obj.Signal, type_: str, figure_folder) -> ca.analysis.peak_picking.ResultPeaks:
    if type_ == "fid":
        signal.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.3), invert=True))
    signal.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True))

    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(signal, scipy_kwargs={"height": 20000, "width": 0.1})
    peak_result_fid = ca.analysis.boundary_detection.rolling_ball(peak_locations, n=10, min_height=0.002,
                                                                  n_points_with_pos_slope=1)

    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # ca.plotting.baseline(signal, fig=fig)
    ca.plotting.signal(signal, fig=fig)
    ca.plotting.peaks(peak_result_fid, fig=fig)

    fig.write_html(figure_folder / f"{type_}_{signal.name}.html")

    return peak_result_fid


def get_compound_data(peak_results, reference_data, concentrations, figure_folder, type_: str = "fid"):
    number_of_compounds = len(reference_data)
    # get peak areas
    areas = np.empty((len(peak_results), number_of_compounds))
    for i, fid_result in enumerate(peak_results):
        areas[i, :] = [peak.stats.area for peak in fid_result.peaks]

    # compute response factors
    response = np.ones(number_of_compounds)
    TCB_index = [chem["label"] == "TCB" for chem in reference_data].index(True)
    areas = (areas.T/areas[:, TCB_index]).T
    for i in range(number_of_compounds):
        if i == TCB_index:
            continue
        x = concentrations[:, i]
        y = areas[:, i]
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        if r_value ** 2 < 0.98:
            print(f"!!! R**2 is large on {reference_data[i]['name']}.  R**2: {r_value ** 2}  !!!")
        response[i] = slope

        fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
        fig.add_trace(go.Scatter(x=x, y=y, mode="markers"))
        x = np.linspace(0, np.max(x))
        fig.add_trace(go.Scatter(x=x, y=slope*x+intercept, mode="lines"))
        fig.layout.title = f"{reference_data[i]['name']}: slope:{slope}, intercept:{intercept},  R**2: {r_value ** 2}"
        fig.write_html(figure_folder / f"{type_}_{reference_data[i]['name']}.html")

    # get retention time
    retention_time = [peak.max_loc for i, peak in enumerate(peak_results[3].peaks)]

    return retention_time, response


def get_compound_data_ms(peak_results, reference_data, concentrations, figure_folder):
    ms_extracted = [ca.analysis.ms_extraction.ms_extract_index(peak.parent, peak.bounds) for peak in
                    peak_results[0].peaks]
    retention_time, response = get_compound_data(peak_results, reference_data, concentrations, figure_folder, "ms")
    return retention_time, response, ms_extracted


def process_group(reference_data, root_folder, pattern, concentrations):
    figure_folder = pathlib.Path(root_folder) / "figs" / pattern.split("-")[2]
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    # get data
    specific_folders = glob.glob(root_folder + pattern)
    specific_folders.sort(key=lambda x: x.split("//")[-1].split("-")[3])
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]

    # process data
    fid_peak_results = []
    ms_peak_results = []
    for ms, fid in data:
        fid_peak_results.append(process_one_signal(fid, "fid", figure_folder))
        ms_peak_results.append(process_one_signal(ms, "ms", figure_folder))

    # check data processed ok
    if any(len(result) != len(reference_data) for result in fid_peak_results):
        raise ValueError(f"Number of peaks detected in fid does not match reference data.")
    if any(len(result) != len(reference_data) for result in ms_peak_results):
        raise ValueError(f"Number of peaks detected in ms does not match reference data.")

    # get processed data
    retention_time_fid, response_fid = get_compound_data(fid_peak_results, reference_data, concentrations, figure_folder)
    retention_time_ms, response_ms, ms_extracted = get_compound_data_ms(ms_peak_results, reference_data, concentrations, figure_folder)

    # add to library
    for i, chem in enumerate(reference_data):
        fid_response = ca.gc_lc.CompoundResponse(
            method="decane_fid",
            retention_time=retention_time_fid[i],
            response=response_fid[i]
        )
        ms_response = ca.gc_lc.CompoundResponse(
            method="decane_ms",
            retention_time=retention_time_ms[i],
            response=response_ms[i],
            mass_spectrum=ms_extracted[i].to_numpy(reduce=True)
        )
        chem["responses"] = [fid_response, ms_response]
        if LIBRARY.find_by_name(chem["name"]):
            LIBRARY.delete_compound(chem["name"])
        LIBRARY.add_compound(
            ca.gc_lc.Compound(
                **chem
            )
        )

    LIBRARY.to_JSON(lib_path, overwrite=True, binary=True)


def run_kn_2():
    reference_data = [
        {
            "label": "K6_2",
            "groups": ["ketone", "methyl_ketone"],
            "name": "2-hexanone",
            "cas": "591-78-6",
            "smiles": "CCCCC(=O)C",
            "density": 0.812,
            "boiling_temperature": 127
        },
        {
            "label": "K7_2",
            "groups": ["ketone", "methyl_ketone"],
            "name": "2-heptanone",
            "cas": "110-43-0",
            "smiles": "CCCCCC(=O)C",
            "density": 0.82,
            "boiling_temperature": 149
        },
        {
            "label": "K8_2",
            "groups": ["ketone", "methyl_ketone"],
            "name": "2-octanone",
            "cas": "111-13-7",
            "smiles": "CCCCCCC(=O)C",
            "density": 0.82,
            "boiling_temperature": 173
        },
        {
            "label": "K9_2",
            "groups": ["ketone", "methyl_ketone"],
            "name": "2-nonanone",
            "cas": "821-55-6",
            "smiles": "CCCCCCCC(=O)C",
            "density": 0.82,
            "boiling_temperature": 192
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        }
    ]
    concentrations = np.array([
        [1.671603295, 1.706956473, 1.82286301, 1.72386927, 1],
        [0.8358016476, 0.8534782363, 0.9114315052, 0.8619346348, 1],
        [0.334320659, 0.3413912945, 0.3645726021, 0.3447738539, 1],
        [0.1671603295, 0.1706956473, 0.182286301, 0.172386927, 1],
    ])
    # concentrations = concentrations.T

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-Kn_2-*.D"

    process_group(reference_data, root_folder, pattern, concentrations)


def main_first():
    lib_ = ca.gc_lc.GCLibrary('DJW-oxidation')
    lib_.to_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON")


if __name__ == "__main__":
    # main_first()
    run_kn_2()
