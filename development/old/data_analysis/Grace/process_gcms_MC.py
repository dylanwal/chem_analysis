import pathlib
from collections import OrderedDict

import numpy as np
import polars as pl

import chem_analysis as ca

lib = ca.library.Library.from_JSON(
        r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json")
rt_lib = ca.a.gc_lc_analysis.RetentionTimeLibrary.from_library(lib, "decane_fid")
lib_df = lib.to_dataframe()


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


def process_one(data_path: pathlib.Path, label: str):
    ms, fid = load_single_file(data_path / "GCMS", label)
    fid = ca.p.baseline.MorphologicalAverage(window=500).run(fid)
    peaks, bounds = ca.a.peaks.find_peaks_and_bounds(
        fid,
        peak_det=ca.a.peaks.PeakDerivative(),
        filters=[
            ca.a.peaks.FilterSpans([5.5, 44]),
            ca.a.peaks.FilterHeight(min_abs=8_000),
            ca.a.peaks.FilterSpacing(min_=0.1),
        ],
        bound_det=ca.a.peaks.BoundFirstIncreaseZero(),#smoother=ca.p.smoothing.Gaussian(std=3)),
    )
    I = ca.a.peaks.integrate_trapz(fid, bounds)
    labels_rt = ca.a.gc_lc.search_by_retention(rt_lib, fid.x[peaks], ca.a.gc_lc.RTMethodNearestN(3, 0.3))

    # ms compare
    # ms_bounds = translate_bounds(bounds, fid.x, ms.x, offset=29.675 - 29.58)
    # ms_peak_data = ca.a.ms_analysis.ms_extract_index(ms, ms_bounds)
    #
    # labels_ms = []
    # for labels_, ms_ in zip(labels_rt, ms_peak_data):
    #     labels_ms.append(
    #         ca.a.ms.search_by_ms_chemicals(
    #             labels_,
    #             ms_,
    #             scorer_ms=ca.a.ms.similarity.earth_movers_distance,
    #             filter_ms=[ca.a.ms_analysis.FilterTopNMatches(n=1)]
    #         )
    #     )

    labels = [l[0] if l else None for l in labels_rt]
    label_dict = OrderedDict()
    for i, label in enumerate(labels):
        if label is None:
            label_dict[str(i)] = label
        else:
            abbr = label.get_identifier(ca.library.identifiers.UserDefined, class_="shorthand").value
            label_dict[abbr] = label
    results = pl.DataFrame({"label": list(label_dict.keys()), "I": I, "bounds": bounds})

    return fid, results, label_dict


def main_baseline():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    times = np.array([0, 30, 60, 120, 240])
    data_label = "DJW-11-53-run44"
    labels = [data_label + f"-t{i}" for i in times]

    figs = []
    for label in labels:
        ms, fid = load_single_file(data_path / "GCMS", label)
        proc = ca.p.Processor(
            ca.p.baseline.MorphologicalAverage(window=500, save_result=True)
        )
        proc_sig = proc.run(fid)
        print(f"window: {proc[0].window}")
        fig = ca.plot.signal(fid)
        fig = ca.plot.baseline(proc, fig=fig)
        fig = ca.plot.signal(proc_sig, fig=fig)
        fig.layout.yaxis.range = [-3_000, 30_000]
        fig.layout.title = "FID"

        figs.append(fig)
        print(f"processed: {label}")

    from chem_analysis.plotting.plotly_plots.plotly_utils import merge_figures
    merge_figures(figs, filename=data_path / 'baseline_fid.html')


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    times = np.array([0, 15, 30, 60, 120, 240])
    data_label = "DJW-11-53-run44"
    labels = [data_label + f"-t{i}" for i in times]

    data = [process_one(data_path, label) for label in labels]

    # re-organizing data
    df = pl.concat([d[1].select("label") for d in data]).unique()
    fid_figs = []
    chemical_dict = dict()
    for t, (fid, df_, chem_dict) in zip(times, data):
        fig = ca.plot.signal(fid)
        fig = ca.plot.add_peaks(fid, df_.get_column("bounds").to_numpy(), mode=[0, 1], fig=fig, labels=df_.get_column("label").to_list())
        fid_figs.append(fig)
        df_sub = df_.rename({"I": f"{t}"})
        df_sub = df_sub.drop("bounds")
        df = df.join(df_sub, on="label", how="left")
        chemical_dict.update(chem_dict)

    # add mM and carbon
    df = df.sort(by="label")
    df = df.fill_null(0)
    rf = lib_df.select(["a_ResponseFactor_c_Method_decane_fid", "i_UserDefined_shorthand"])
    rf = rf.rename({"i_UserDefined_shorthand": "label", "a_ResponseFactor_c_Method_decane_fid": "RF"})
    df = df.join(rf, on="label", how="left")
    carbon = pl.DataFrame(
        {
            "label": list(chemical_dict.keys()),
            "carbon": list(c.get_identifier("ChemicalFormula").element_count("C") - 3*c.get_identifier("ChemicalFormula").element_count("Si") for c in chemical_dict.values())
        }
    )
    df = df.join(carbon, on="label", how="left")
    TCB_row = df[-1]
    mmol_TCB = 0.055  # M
    decane_conc = 0.777 / 142.1 * 1000 / 15  # (.777 g / 142.27 g/mol * 1000 to make mmol) / 15 ml of rxn vol * 0.025 ml = moles of decane in vial
    df = df.with_columns((pl.col(col)/TCB_row[col]*mmol_TCB*pl.col("RF")).alias(f"{col}_conc") for col in df.columns[1:-2])
    df = df.with_columns((pl.col(col)*pl.col("carbon")).alias(f"{col}_C") for col in df.columns[1:-1] if col.endswith("_conc"))

    # saving data
    df.write_csv(data_path / (data_label + '_data.csv'))
    pl.Config.set_tbl_rows(100)
    print(df)

    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=data_path / (data_label + '_fid.html'))


def main_plot():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    times = np.array([0, 30, 60, 120, 240])
    data_label = "DJW-11-53-run44"
    labels = [data_label + f"-t{i}" for i in times]

    sigs = [load_single_file(data_path / "GCMS", label) for label in labels]
    fids = ca.gc_lc.GCSignal2D.from_signals([i[1] for i in sigs], times, "time",
                                            unify_method=ca.base_obj.UnifyMethodExpandInterpolate())
    mss = [i[0] for i in sigs]

    fig = ca.plot.signal2D_overlap_signals(fids)
    fig.show()


if __name__ == "__main__":
    # main_baseline()
    main()
    # main_plot()
