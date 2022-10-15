import glob
import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)

wells = [
    "A1",
    "A2",
    "A9",
    "A16",
    "A23",
    "A24",
    "B1",
    "B2",
    "B23",
    "B24",
    "C11",
    "C14",
    "D4",
    "D21",
    "E12",
    "E13",
    "F6",
    "F19",
    "G10",
    "G15",
    "H1",
    "H8",
    "H12",
    "H13",
    "H17",
    "H24",
    "I1",
    "I8",
    "I12",
    "I13",
    "I17",
    "I24",
    "J10",
    "J15",
    "K6",
    "K19",
    "L12",
    "L13",
    "M4",
    "M21",
    "N11",
    "N14",
    "O1",
    "O2",
    "O23",
    "O24",
    "P1",
    "P2",
    "P9",
    "P16",
    "P23",
    "P24",

]


def main():
    folder = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Instrument\polymerization\RAFT4-10min"
    files = glob.glob(folder + "\\*.csv")
    data = dict()
    for file in files:
        key = file.split("\\")[-1].replace(".csv", "")
        with open(file, 'r') as f:
            data[key] = pd.read_csv(f)

    def cal_func_RI(time: np.ndarray):
        return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    signals = []
    for name, df in data.items():
        signals.append(chem.SECSignal(name=name, xy=df[["time_RI", "RI"]].to_numpy(), cal=cal_RI,
                                      x_label="retention time (min)", y_label="signal"))

    count = 0
    df = None
    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline)
        sig.peak_picking(lb=10.85, ub=12.8)   # lb=10.85, ub=12.2

        if count == 0:
            count += 1
            df = sig.stats_df()
            sig.plot()
        else:
            df = pd.concat([df, sig.stats_df()], axis=1)

    df = df.T
    df = df[:52]
    df.index = wells

    import well_plate
    wp = well_plate.WellPlate(384, "rect")
    wp.add_data(df["mw_n"])
    wp.plot(key="mw_n")
    print("mean :", np.mean(df["mw_n"]))
    print("std: ", np.std(df["mw_n"]))

    print("hi")


if __name__ == '__main__':
    main()
