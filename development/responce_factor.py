import glob
import pathlib
import shutil

import numpy as np
from scipy import stats
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)


def process_one_signal(signal: ca.base_obj.Signal, type_: str, figure_folder,
                       peak_mask) -> ca.analysis.peak_picking.ResultPeaks:
    if type_ == "fid":
        signal.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.7), invert=True))
    signal.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True))

    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(signal,
                                                               mask=peak_mask,
                                                               scipy_kwargs={"height": 6000, "width": 0.1}
                                                               )
    peak_result_fid = ca.analysis.boundary_detection.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                                  n_points_with_pos_slope=2)

    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # ca.plotting.baseline(signal, fig=fig)
    ca.plotting.signal(signal, fig=fig)
    ca.plotting.peaks(peak_result_fid, fig=fig)

    fig.write_html(figure_folder / f"signal_{type_}_{signal.name}.html", include_plotlyjs='cdn')

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
    areas = (areas.T / areas[:, TCB_index]).T
    for i in range(number_of_compounds):
        if i == TCB_index:
            continue
        x = areas[:, i]
        y = concentrations[:, i]
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        if r_value ** 2 < 0.98:
            print(f"!!! R**2 is large on {reference_data[i]['name']}.  R**2: {r_value ** 2}  !!!")
        response[i] = slope

        fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
        fig.add_scatter(x=x, y=y, mode="markers")
        x = np.linspace(0, np.max(x))
        fig.add_scatter(x=x, y=slope * x + intercept, mode="lines")
        fig.layout.title = f"{reference_data[i]['name']}: slope:{slope}, intercept:{intercept},  R**2: {r_value ** 2}"
        fig.write_html(figure_folder / f"fit_{type_}_{reference_data[i]['name']}.html", include_plotlyjs='cdn')

    # get retention time
    retention_time = [peak.max_loc for i, peak in enumerate(peak_results[3].peaks)]

    return retention_time, response


def get_compound_data_ms(peak_results, reference_data, concentrations, figure_folder):
    ms_extracted = [ca.analysis.ms_extraction.ms_extract_index(peak.parent, peak.bounds) for peak in
                    peak_results[0].peaks]
    retention_time, response = get_compound_data(peak_results, reference_data, concentrations, figure_folder, "ms")

    for i, ms in enumerate(ms_extracted):
        fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
        ca.plotting.signal(ms, fig=fig)
        fig.layout.title = f"{reference_data[i]['name']}"
        fig.write_html(figure_folder / f"ms_{reference_data[i]['name']}.html", include_plotlyjs='cdn')
    return retention_time, response, ms_extracted


def process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask):
    figure_folder = pathlib.Path(root_folder) / "figs" / figure_label
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    # get data
    specific_folders = glob.glob(root_folder + pattern)
    if len(specific_folders) != 4:
        raise ValueError("Not correct number of files found.")
    specific_folders.sort(key=lambda x: x.split("//")[-1].replace("redo", "").replace("-", "").replace("_", ""))
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]

    # process data
    fid_peak_results = []
    ms_peak_results = []
    for ms, fid in data:
        fid_peak_results.append(process_one_signal(fid, "fid", figure_folder, peak_mask))
        ms_peak_results.append(process_one_signal(ms, "ms", figure_folder, peak_mask))

    # check data processed ok
    if any(len(result) != len(reference_data) for result in fid_peak_results) or not fid_peak_results:
        raise ValueError(f"Number of peaks detected in fid does not match reference data.")
    if any(len(result) != len(reference_data) for result in ms_peak_results) or not ms_peak_results:
        raise ValueError(f"Number of peaks detected in ms does not match reference data.")

    # get processed data
    retention_time_fid, response_fid = get_compound_data(fid_peak_results, reference_data, concentrations,
                                                         figure_folder)
    retention_time_ms, response_ms, ms_extracted = get_compound_data_ms(ms_peak_results, reference_data, concentrations,
                                                                        figure_folder)

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


def run_kn_2():
    reference_data = [  # order of chemical must match chromatogram
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

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-Kn_2-[0-9].D"
    figure_label = "Kn_2"
    peak_mask = ca.processing.weigths.Spans([4, 31], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_k10():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "K10_3",
            "groups": ["ketone"],
            "name": "3-decanone",
            "cas": "928-80-3",
            "smiles": "CCCCCCCC(=O)CC",
            "density": 0.825,
            "boiling_temperature": 204
        },
        {
            "label": "K10_2",
            "groups": ["ketone"],
            "name": "2-decanone",
            "cas": "693-54-9",
            "smiles": "CCCCCCCC(=O)CC",
            "density": 0.825,
            "boiling_temperature": 211
        },
        {
            "label": "K10",
            "groups": ["ketone"],
            "name": "decanal",
            "cas": "112-31-2",
            "smiles": "CCCCCCCC(=O)CC",
            "density": 0.825,
            "boiling_temperature": 207
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
    ]
    concentrations = np.array([
        [1.693229322, 1.684864746, 1.732662327, 1],
        [0.8466146612, 0.8424323728, 0.8663311635, 1],
        [0.3386458645, 0.3369729491, 0.3465324654, 1],
        [0.1693229322, 0.1684864746, 0.1732662327, 1]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-K_[0-9].D"
    figure_label = "K10"
    peak_mask = ca.processing.weigths.Spans([[24.5, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_k10_4():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "K10_4",
            "groups": ["ketone"],
            "name": "4-decanone",
            "cas": "624-16-8",
            "smiles": "CCCCCCC(=O)CCC",
            "density": 0.825,
            "boiling_temperature": 208
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
    ]

    concentrations = np.array([
        [2.274374526, 1],
        [0.9097498103, 1],
        [0.4548749051, 1],
        [0.1819499621, 1]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-K10_4_[0-9]_solo.D"
    figure_label = "K10_4"
    peak_mask = ca.processing.weigths.Spans([[20, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_k10_5():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "K10_5",
            "groups": ["ketone"],
            "name": "5-decanone",
            "cas": "820-29-1",
            "smiles": "CCCCCC(=O)CCCC",
            "density": 0.825,
            "boiling_temperature": 210
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
    ]
    concentrations = np.array([
        [2.430199608, 1],
        [0.9720798431, 1],
        [0.4860399215, 1],
        [0.1944159686, 1]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-K10_5_[0-9]_solo.D"
    figure_label = "K10_5"
    peak_mask = ca.processing.weigths.Spans([[20, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_da():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
        {
            "label": "DA4",
            "groups": ["dicarboxylic acid", "TMS"],
            "name": "succinic acid, 2 x TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCC(=O)O[Si](C)(C)C",
            "parent": "succinic acid"
        },
        {
            "label": "DA5",
            "groups": ["dicarboxylic acid", "TMS"],
            "name": "glutaric acid, 2 x TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCC(=O)O[Si](C)(C)C",
            "parent": "glutaric acid"
        },
        {
            "label": "DA6",
            "groups": ["dicarboxylic acid", "TMS"],
            "name": "adipic acid, 2 x TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCC(=O)O[Si](C)(C)C",
            "parent": "adipic acid"
        },
        {
            "label": "DA7",
            "groups": ["dicarboxylic acid", "TMS"],
            "name": "pimelic acid, 2 x TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCC(=O)O[Si](C)(C)C",
            "parent": "pimelic acid"
        },
        {
            "label": "DA8",
            "groups": ["dicarboxylic acid", "TMS"],
            "name": "suberic acid, 2 x TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCCC(=O)O[Si](C)(C)C",
            "parent": "suberic acid"
        }
    ]
    concentrations = np.array([
        [1, 1.741610267, 1.707593705, 1.549707538, 1.777826902, 1.828741261],
        [1, 0.8708051335, 0.8537968527, 0.7748537692, 0.888913451, 0.9143706305],
        [1, 0.3483220534, 0.3415187411, 0.3099415077, 0.3555653804, 0.3657482522],
        [1, 0.1741610267, 0.1707593705, 0.1549707538, 0.1777826902, 0.1828741261],

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-DA[0-9]-TMS.D"
    figure_label = "DA"
    peak_mask = ca.processing.weigths.Spans([[28, 43.03], [43.3, 44.8]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_ca1():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "CA3",
            "groups": ["carboxylic acid", "TMS"],
            "name": "propionic acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CC",
            "parent": "propionic acid"
        },
        {
            "label": "CA4",
            "groups": ["carboxylic acid", "TMS"],
            "name": "butyric acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCC",
            "parent": "butyric acid"
        },
        {
            "label": "CA6",
            "groups": ["carboxylic acid", "TMS"],
            "name": "hexanoic acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCC",
            "parent": "hexanoic acid"
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
        {
            "label": "CA10",
            "groups": ["carboxylic acid", "TMS"],
            "name": "decanoic acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCCCCC",
            "parent": "decanoic acid"
        }
    ]
    concentrations = np.array([
        [1.741494412, 1.741494412, 1.733265358, 1, 1.741494412],
        [0.6965977648, 0.6965977648, 0.6933061431, 1, 0.6965977648],
        [0.1741494412, 0.1741494412, 0.1733265358, 1, 0.1741494412],
        [0.06965977648, 0.06965977648, 0.06933061431, 1, 0.06965977648]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-CA1_[0-9].D"
    figure_label = "CA"
    peak_mask = ca.processing.weigths.Spans([[None, 4], [5.5, 7], [14, 18], [28, 31], [38, 41]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_ca2():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "CA5",
            "groups": ["carboxylic acid", "TMS"],
            "name": "valeric acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCC",
            "parent": "valeric acid"
        },
        {
            "label": "CA7",
            "groups": ["carboxylic acid", "TMS"],
            "name": "heptanic acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCC",
            "parent": "heptanic acid"
        },
        {
            "label": "CA8",
            "groups": ["carboxylic acid", "TMS"],
            "name": "octanoic acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCCC",
            "parent": "octanoic acid"
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
        {
            "label": "CA9",
            "groups": ["carboxylic acid", "TMS"],
            "name": "nonanioc acid, TMS",
            "cas": "",
            "smiles": "[Si](C)(C)(C)OC(=O)CCCCCCCC",
            "parent": "nonanioc acid"
        }
    ]
    concentrations = np.array([
        [1.71140894, 1.733734738, 1.730610167, 1, 1.745115273],
        [0.8557044701, 0.8668673688, 0.8653050837, 1, 0.8725576363],
        [0.3422817881, 0.3467469475, 0.3461220335, 1, 0.3490230545],
        [0.171140894, 0.1733734738, 0.1730610167, 1, 0.1745115273]
    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-CA2_[0-9].D"
    figure_label = "CA2"
    peak_mask = ca.processing.weigths.Spans([[8, 27.8], [29, 35]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_C10():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "C10",
            "groups": ["alkane"],
            "name": "decane",
            "cas": "124-18-5",
            "smiles": "CCCCCCCCCC",
            "density": 0.73,
            "boiling_temperature": 174
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
    ]

    concentrations = np.array([
        [3.433293764, 1],
        [1.716646882, 1],
        [0.6866587528, 1],
        [0.3433293764, 1]
    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal_D_[0-9].D"
    figure_label = "C10"
    peak_mask = ca.processing.weigths.Spans([[5, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_A10():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "A10_5",
            "groups": ["alcohol", "TMS"],
            "name": "5-decanol, TMS",
            "cas": "",
            "smiles": "CCCCCC(O[Si](C)(C)C)CCCC",
            "parent": "5-decanol"
        },
        {
            "label": "A10_4",
            "groups": ["alcohol", "TMS"],
            "name": "4-decanol, TMS",
            "cas": "",
            "smiles": "CCCCCCC(O[Si](C)(C)C)CCC",
            "parent": "4-decanol"
        },
        {
            "label": "A10_3",
            "groups": ["alcohol", "TMS"],
            "name": "3-decanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCC(O[Si](C)(C)C)CC",
            "parent": "3-decanol"
        },
        {
            "label": "A10_2",
            "groups": ["alcohol", "TMS"],
            "name": "2-decanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCCC(O[Si](C)(C)C)C",
            "parent": "2-decanol"
        },
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
        {
            "label": "A10",
            "groups": ["alcohol", "TMS"],
            "name": "1-decanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCCCC(O[Si](C)(C)C)",
            "parent": "1-decanol"
        },
    ]
    concentrations = np.array([
        [0.3837045671, 1.711547729, 1.697907536, 1.690197861, 1, 1.702651951],
        [0.1918522835, 0.8557738644, 0.8489537678, 0.8450989305, 1, 0.8513259753],
        [0.07674091341, 0.3423095458, 0.3395815071, 0.3380395722, 1, 0.3405303901],
        [0.03837045671, 0.1711547729, 0.1697907536, 0.1690197861, 1, 0.1702651951]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-A10-TMS-[0-9].D"
    figure_label = "A10"
    peak_mask = ca.processing.weigths.Spans([[24.3, 25.8], [26.1, 33]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_A1():
    reference_data = [  # order of chemical must match chromatogram
        # {
        #     "label": "A3_1",
        #     "groups": ["alcohol", "TMS"],
        #     "name": "1-propanol, TMS",
        #     "cas": "",
        #     "smiles": "CCCO[Si](C)(C)C",
        #     "parent": "1-propanol"
        # },
        # {
        #     "label": "A4_1",
        #     "groups": ["alcohol", "TMS"],
        #     "name": "1-butanol, TMS",
        #     "cas": "",
        #     "smiles": "CCCCO[Si](C)(C)C",
        #     "parent": "1-butanol"
        # },
        {
            "label": "A5_1",
            "groups": ["alcohol", "TMS"],
            "name": "1-pentanol, TMS",
            "cas": "",
            "smiles": "CCCCCO[Si](C)(C)C",
            "parent": "1-pentanol"
        },
        {
            "label": "A6_1",
            "groups": ["alcohol", "TMS"],
            "name": "1-hexanol, TMS",
            "cas": "",
            "smiles": "CCCCCCO[Si](C)(C)C",
            "parent": "1-hexanol"
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
        [1.748517449, 1.808756319, 1.751107878, 1.748612949, 1],
        [0.8742587244, 0.9043781595, 0.8755539389, 0.8743064746, 1],
        [0.3497034898, 0.3617512638, 0.3502215755, 0.3497225898, 1],
        [0.1748517449, 0.1808756319, 0.1751107878, 0.1748612949, 1]
    ])
    concentrations = concentrations[:, 2:]

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-A1-[0-9]-TMS.D"
    figure_label = "A1"
    peak_mask = ca.processing.weigths.Spans([[5, 9.5], [28, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_A2():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "A7_1",
            "groups": ["alcohol", "TMS"],
            "name": "1-heptanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCO[Si](C)(C)C",
            "parent": "1-heptanol"
        },
        {
            "label": "A8_1",
            "groups": ["alcohol", "TMS"],
            "name": "1-octanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCCO[Si](C)(C)C",
            "parent": "1-octanol"
        },
        {
            "label": "A9_1",
            "groups": ["alcohol", "TMS"],
            "name": "1-nonanol, TMS",
            "cas": "",
            "smiles": "CCCCCCCCCO[Si](C)(C)C",
            "parent": "1-nonanol"
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
        [1.904742484, 1.763449104, 1.732274244, 1],
        [0.9523712422, 0.881724552, 0.866137122, 1],
        [0.3809484969, 0.3526898208, 0.3464548488, 1],
        [0.1904742484, 0.1763449104, 0.1732274244, 1]
    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-A2-[0-9]-TMS.D"
    figure_label = "A2"
    peak_mask = ca.processing.weigths.Spans([[13.5, 14.5], [19, 21], [25, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_AA10():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "TCB",
            "groups": "standard",
            "name": "1,2,3-trichlorobenzene",
            "cas": "87-61-6",
            "smiles": "C1=CC(=C(C(=C1)Cl)Cl)Cl",
            "density": 1.45,
            "boiling_temperature": 218.5
        },
        {
            "label": "AA10_5",
            "groups": ["alcohol", "acetylated"],
            "name": "5-decanol, acetylated",
            "cas": "",
            "smiles": "CCCCCC(OC(=O)C)CCCC",
            "parent": "5-decanol"
        },
        {
            "label": "AA10_4",
            "groups": ["alcohol", "acetylated"],
            "name": "4-decanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCC(OC(=O)C)CCC",
            "parent": "4-decanol"
        },
        {
            "label": "AA10_3",
            "groups": ["alcohol", "acetylated"],
            "name": "3-decanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCC(OC(=O)C)CC",
            "parent": "3-decanol"
        },
        {
            "label": "AA10_2",
            "groups": ["alcohol", "acetylated"],
            "name": "2-decanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCC(OC(=O)C)C",
            "parent": "2-decanol"
        },
        {
            "label": "AA10",
            "groups": ["alcohol", "acetylated"],
            "name": "1-decanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCCCCOC(=O)C",
            "parent": "1-decanol"
        },
    ]
    concentrations = np.array([
        [1, 0.3837045671, 1.711547729, 1.697907536, 1.690197861, 1.702651951],
        [1, 0.1918522835, 0.8557738644, 0.8489537678, 0.8450989305, 0.8513259753],
        [1, 0.07674091341, 0.3423095458, 0.3395815071, 0.3380395722, 0.3405303901],
        [1, 0.03837045671, 0.1711547729, 0.1697907536, 0.1690197861, 0.1702651951]

    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-AA10-[0-9]-redo.D"
    figure_label = "AA10"
    peak_mask = ca.processing.weigths.Spans([[28, 32.6], [33.45, 40]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_AA1():
    reference_data = [  # order of chemical must match chromatogram
        # {
        #     "label": "AA3_1",
        #     "groups": ["alcohol", "acetylated"],
        #     "name": "1-propanol, acetylated",
        #     "cas": "",
        #     "smiles": "CCCOC(=O)C",
        #     "parent": "1-propanol"
        # },
        {
            "label": "AA4_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-butanol, acetylated",
            "cas": "",
            "smiles": "CCCCOC(=O)C",
            "parent": "1-butanol"
        },
        {
            "label": "AA5_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-pentanol, acetylated",
            "cas": "",
            "smiles": "CCCCCOC(=O)C",
            "parent": "1-pentanol"
        },
        {
            "label": "AA6_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-hexanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCOC(=O)C",
            "parent": "1-hexanol"
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
        [1.748517449, 1.808756319, 1.751107878, 1.748612949, 1],
        [0.8742587244, 0.9043781595, 0.8755539389, 0.8743064746, 1],
        [0.3497034898, 0.3617512638, 0.3502215755, 0.3497225898, 1],
        [0.1748517449, 0.1808756319, 0.1751107878, 0.1748612949, 1]
    ])
    concentrations = concentrations[:, 1:]

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-AA1-[0-9].D"
    figure_label = "AA1"
    peak_mask = ca.processing.weigths.Spans([[5, 6], [8, 9.4], [12, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def run_AA2():
    reference_data = [  # order of chemical must match chromatogram
        {
            "label": "A7_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-heptanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCOC(=O)C",
            "parent": "1-heptanol"
        },
        {
            "label": "A8_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-octanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCCOC(=O)C",
            "parent": "1-octanol"
        },
        {
            "label": "A9_1",
            "groups": ["alcohol", "acetylated"],
            "name": "1-nonanol, acetylated",
            "cas": "",
            "smiles": "CCCCCCCCCOC(=O)C",
            "parent": "1-nonanol"
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
        [1.904742484, 1.763449104, 1.732274244, 1],
        [0.9523712422, 0.881724552, 0.866137122, 1],
        [0.3809484969, 0.3526898208, 0.3464548488, 1],
        [0.1904742484, 0.1763449104, 0.1732274244, 1]
    ])

    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-A2-[0-9].D"
    figure_label = "AA2"
    peak_mask = ca.processing.weigths.Spans([[10, 24.5], [28, 31]], invert=True)

    process_group(reference_data, root_folder, pattern, concentrations, figure_label, peak_mask)


def main_first():
    lib_ = ca.gc_lc.GCLibrary('DJW-oxidation')
    lib_.to_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON")


if __name__ == "__main__":
    # main_first()
    # run_kn_2()
    # run_k10()
    # run_k10_4()
    # run_k10_5()
    # run_da()
    # run_ca1()
    # run_ca2()
    # run_C10()
    # run_A10()
    # run_A1()
    # run_A2()
    # run_AA10()
    # run_AA1()
    # run_AA2()  # not the best data

    run_k_1()
    run_HA()
    run_Br10()
    run_lactone()
    run_BHA() # branched hydroxy acids

    # LIBRARY.to_JSON(lib_path, overwrite=True, binary=True)
