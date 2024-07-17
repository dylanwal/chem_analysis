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


def main():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\10\10_10")
    fid_timeseries = ca.gc_lc.GCSignal2D.from_npz(data_path / "FID_timeseries.npz")
    ms_timeseries = ca.gc_lc.GCMSSignal2D.from_npz(data_path / "MS_timeseries.npz")
    print(fid_timeseries)
    print(ms_timeseries)


if __name__ == "__main__":
    main()
