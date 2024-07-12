import pathlib
import re
import shutil

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca
from development.time_series_support import ResultTimeSeries, plot_results

# lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
# LIBRARY = ca.gc_lc.GCLibrary.from_JSON(lib_path)
# INTERNAL_STANDARD = LIBRARY.find_by_label("TCB")
# picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
# picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")
# PLOTTING_GROUPS = ["dicarboxylic acid", "hydroxy acids", "carboxylic acid", "alcohol", "methyl_ketone", "ketone", "peroxide", "alkane"]


def get_folders(data_path: pathlib.Path, pattern: str) -> tuple[list[pathlib.Path], np.ndarray]:
    pattern_compile = re.compile(pattern.replace("*", "([0-9]+)"))

    def sort_func(x: pathlib.Path) -> int:
        return int(pattern_compile.match(x.name).groups()[0])

    specific_folders = list(data_path.glob(pattern))
    print(len(specific_folders), "files found for analysis")

    specific_folders.sort(key=sort_func)
    times_ = [sort_func(folder) for folder in specific_folders]

    return [data_path / file for file in specific_folders], np.array(times_)


def process_timeseries(data_path: str, pattern: str):
    if isinstance(data_path, str):
        data_path = pathlib.Path(data_path)
    folders, times = get_folders(data_path, pattern)

    # get data
    fids, mss = [], []
    for folder in folders:
        ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(folder)

    fids = ca.gc_lc.GCMSSignal2D.from_signals(fids, times)
    mss = ca.gc_lc.GCMSSignal2D.from_signals(mss, times)


    # fid
    fids.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 4), invert=True))
    fids.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )

    # ms
    mss.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )

    ca.plotting.qt(fids)


def main():
    data_path = r"C:\Users\nicep\Desktop\research_wis\data\11\11_29\GC_MS"
    pattern = "DJW-11-29-*min.D"
    process_timeseries(data_path, pattern)


if __name__ == "__main__":
    main()
