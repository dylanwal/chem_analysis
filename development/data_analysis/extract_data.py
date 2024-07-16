import pathlib
import re

import numpy as np

import chem_analysis as ca


def get_folders(data_path: str | pathlib.Path, pattern: str) -> tuple[list[pathlib.Path], np.ndarray]:
    if isinstance(data_path, str):
        data_path = pathlib.Path(data_path)

    pattern_compile = re.compile(pattern.replace("*", "([0-9]+)"))

    def sort_func(x: pathlib.Path) -> int:
        result = pattern_compile.match(x.name)
        if result is None:
            raise ValueError(f"'{x.name}' is not a valid pattern for '{pattern}'.")
        return int(result.groups()[0])

    specific_folders = list(data_path.glob(pattern))
    print(len(specific_folders), "files found for analysis")

    specific_folders.sort(key=sort_func)
    times_ = [sort_func(folder) for folder in specific_folders]

    if len(specific_folders) == 0:
        raise RuntimeError("No folders found for analysis")

    return [data_path / file for file in specific_folders], np.array(times_)


def build_timeseries(folders: list[pathlib.Path], times: np.ndarray) \
        -> tuple[ca.gc_lc.GCMSSignal2D, ca.gc_lc.GCMSSignal2D]:
    fids, mss = [], []
    for folder in folders:
        ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(folder)
        fids.append(fid)
        mss.append(ms)
    print("files parsing complete")

    # make fids all the same length
    min_fid_length = min([len(fid.x_raw) for fid in fids])
    for fid in fids:
        fid.x_raw = fid.x_raw[:min_fid_length]
        fid.y_raw = fid.y_raw[:min_fid_length]

    fid_timeseries = ca.gc_lc.GCSignal2D.from_signals(fids, times)
    ms_timeseries = ca.gc_lc.GCMSSignal2D.from_signals(mss, times)
    return fid_timeseries, ms_timeseries


def main():
    data_path = r"C:\Users\nicep\Desktop\10_10"
    folders, times = get_folders(data_path, "DJW-10-10-*h-TMS.D")
    fid, ms = build_timeseries(folders, times)

    fid.to_npz(r"C:\Users\nicep\Desktop\FID_timeseries.npz")
    ms.to_npz(r"C:\Users\nicep\Desktop\MS_timeseries.npz")


if __name__ == "__main__":
    main()
