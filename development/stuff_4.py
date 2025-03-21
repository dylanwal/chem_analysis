import pathlib

import numpy as np
import matplotlib.pyplot as plt

import chem_analysis as ca


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


def translate_bounds(old_bounds: np.ndarray, x_old: np.ndarray, x_new: np.ndarray, offset: float = 0) -> np.ndarray:
    x_left = x_old[old_bounds[:, 0]] + offset
    x_right = x_old[old_bounds[:, 1]] + offset

    left_indices = np.searchsorted(x_new, x_left, side='left')
    right_indices = np.searchsorted(x_new, x_right, side='right') - 1

    return np.column_stack((left_indices, right_indices))


def offset_scorer(mz: np.ndarray, lib: np.ndarray, vec: np.ndarray) -> np.ndarray:
    # shift ms over 1 above the mz 150
    index = np.argmin(np.abs(mz - 150))
    vec_ = np.copy(vec)
    vec_[index - 1:-1] = vec[index:]
    vec_[-1] = 0

    return ca.a.ms_analysis.ScorerDot()(mz, lib, vec_)


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_79")
    data_label = "DJW-11-79-6hr-distilate"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(data_path / (data_label + ".D"))

    fid = ca.p.baseline.MorphologicalAverage().run(fid)
    peaks, bounds = ca.a.peaks.find_peaks_and_bounds(
        fid,
        peak=ca.a.peaks.PeakDerivative(),
        filters=[
            ca.a.peaks.FilterSpans([4, 44]),
            ca.a.peaks.FilterHeight(min_abs=30_000),
            # ca.a.peaks.FilterSpacing(min_=0.3),
        ],
        bounds=ca.a.peaks.BoundFirstIncrease(smoother=ca.p.smoothing.Gaussian(std=5)),
        # (window=40, mode="first_increase_zero"),
        # smoother=ca.p.smoothing.Gaussian(std=5)
    )
    I = ca.a.peaks.integrate_trapz(fid, bounds)

    lib = ca.library.Library.from_JSON(
        r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json")
    # retention time
    rt_lib = ca.a.gc_lc_analysis.RetentionTimeLibrary.from_library(lib, "decane_fid")
    labels_rt = ca.a.gc_lc.search_by_retention(rt_lib, fid.x[peaks], ca.a.gc_lc.RTMethodNearestN(3, 0.3))

    # ms compare
    ms_bounds = translate_bounds(bounds, fid.x, ms.x, offset=29.675 - 29.58)
    ms_peak_data = ca.a.ms_analysis.ms_extract_index(ms, ms_bounds)

    labels_ms = []
    for labels_, ms_ in zip(labels_rt, ms_peak_data):
        labels_ms.append(
            ca.a.ms.search_by_ms_chemicals(
                labels_,
                ms_,
                scorer_ms=ca.a.ms_analysis.ScorerMultiple([ca.a.ms_analysis.ScorerDot(), offset_scorer]),
                filter_ms=[ca.a.ms_analysis.FilterMinScore(0.6), ca.a.ms_analysis.FilterTopNMatches(n=1)]
            )
        )

    # plot
    plt.figure(figsize=(10, 5))

    plt.plot(fid.x, fid.y_normalized_by_max((5,45)))

    plt.xlim(5, 45)
    plt.ylim(-0.1, 1.1)
    # plt.legend(loc="upper right")
    plt.xlabel("Time (min)")
    plt.ylabel("Normalized Intensity")
    plt.show()


if __name__ == "__main__":
    main()
