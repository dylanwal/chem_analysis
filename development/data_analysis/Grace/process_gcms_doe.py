import pathlib
from collections import OrderedDict

import numpy as np
import polars as pl
import plotly.graph_objs as go

import chem_analysis as ca
from chem_analysis.utils.math import get_slice

import development.data_analysis.utils.plotly_utils  # don't delete; automatically adds formating

lib = ca.library.Library.from_JSON(
    r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json")
rt_lib = ca.a.gc_lc_analysis.RetentionTimeLibrary.from_library(lib, "decane_fid")
lib_df = lib.to_dataframe()

INTEGRATE = (
    ('CA4', 5.87, 6.29),
    ('C10', 8.87, 9.406),
    ('CA5', 9.85, 10.41),
    # ('L4_4', 14.91, 15.13),
    ('CA6', 15.13, 15.45),
    # ('L5_4', 16.5, 17.),
    # ('HAA2', 17.67, 18.07),
    ('CA7', 21.22, 21.63),
    # ('L6_4', 23.06, 23.37),
    ('K10_4', 24.04, 24.27),
    ('A10_5', 24.27, 24.49),
    ('A10_4', 24.77, 25.1),
    ('K10_3', 25.32, 25.55),
    ('K10_2', 25.89, 26.12),
    ('A10_3', 26.12, 26.3),
    ('A10_2', 27.02, 27.28),
    ('CA8', 27.46, 27.75),
    ('TCB', 29.41, 30.),
    ('DA4', 32.44, 32.81),
    ('DA5', 38.15, 38.47),
    ('DA6', 42.76, 42.88)
)


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / "GCMS" / f"{label}_data.ms"
    fid_file = data_path / "GCMS" / f"{label}_FID1A.ch"
    ini_file = data_path / "GCMS" / f"{label}_pre_post.ini"
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
    ms, fid = load_single_file(data_path, label)
    fid = ca.p.baseline.MorphologicalAverage(window=500).run(fid)

    bounds = [get_slice(fid.x, start, end) for i, start, end in INTEGRATE]
    bounds = np.array([(s.start, s.stop) for s in bounds])
    I = ca.a.peaks.integrate_trapz(fid, bounds)

    mask = I > 2000
    I = I[mask]
    bounds = bounds[mask]
    labels = tuple(l for i, (l, start, end) in enumerate(INTEGRATE) if mask[i])
    label_dict = {label: lib.find_by_identifier(label, ca.library.i.UserDefined, "shorthand")[0] for label in labels}

    results = pl.DataFrame({"label": labels, "I": I, "bounds": bounds})

    return fid, results, label_dict


def process_one_series(data_path: pathlib.Path, data_label: str, labels: list[str], times: np.ndarray):
    data = [process_one(data_path, label) for label in labels]

    # re-organizing data
    df = pl.concat([d[1].select("label") for d in data]).unique()
    fid_figs = []
    chemical_dict = dict()
    for t, (fid, df_, chem_dict) in zip(times, data):
        fig = ca.plot.signal(fid)
        fig = ca.plot.add_peaks(fid, df_.get_column("bounds").to_numpy(), mode=[0, 1], fig=fig,
                                labels=df_.get_column("label").to_list())
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
            "carbon": list(c.get_identifier("ChemicalFormula").element_count("C") - 3 * c.get_identifier(
                "ChemicalFormula").element_count("Si") for c in chemical_dict.values())
        }
    )
    df = df.join(carbon, on="label", how="left")
    TCB_row = df[-1]
    mmol_TCB = 0.055  # M
    decane_conc = 0.777 / 142.1 * 1000 / 15  # (.777 g / 142.27 g/mol * 1000 to make mmol) / 15 ml of rxn vol * 0.025 ml = moles of decane in vial
    df = df.with_columns(
        (pl.col(col) / TCB_row[col] * mmol_TCB * pl.col("RF")).alias(f"{col}_conc") for col in df.columns if col.isdigit())
    df = df.with_columns(
        (pl.col(col) * pl.col("carbon")).alias(f"{col}_C") for col in df.columns if col.endswith("_conc"))

    # saving data
    df.write_csv(data_path / (data_label + '_data.csv'))
    pl.Config.set_tbl_rows(100)
    print(df)

    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=data_path / (data_label + '_fid.html'))


def main():
    import os
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    files = os.listdir(data_path / "GCMS")
    files = list(filter(lambda x: x.endswith('ini'), files))
    files = list(map(lambda x: x.replace("_pre_post.ini", ""), files))
    chems = set(map(lambda x: x.split("-")[3], files))
    print(files)
    print(chems)

    data = []
    for chem in chems:
        data_label = f"DJW-11-53-{chem}"
        times = np.array([int(f.split("-")[-1][1:]) for f in files if f"-{chem}-" in f])
        times = np.sort(times)
        labels = [data_label + f"-t{t}" for t in times]
        data.append([data_label, labels, times])

    for data_label, labels, times in data:
        process_one_series(data_path, data_label, labels, times)


if __name__ == "__main__":
    main()
