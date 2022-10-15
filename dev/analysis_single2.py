import logging
import functools

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

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
    # load data
    file_name = r"C:\Users\nicep\Downloads\RAFT4-10min.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    for col in list(df.columns):
        if str(col).startswith("Unnamed:"):
            del df[col]

    def cal_func_RI(time: np.ndarray):
            return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    signals = []
    for col in df.columns:
        signals.append(chem.SECSignal(name=col, ser=df[col], cal=cal_RI,
                                      x_label="retention time (min)", y_label="signal"))

    count = 0
    df = None
    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7)
        sig.peak_picking(lb=10.1, ub=12)   # lb=10.85, ub=12.2

        if count == 0:
            count += 1
            df = sig.stats_df()
            sig.plot()
        else:
            df = pd.concat([df, sig.stats_df()], axis=1)

    df = df.T
    # df_ = df[:52]
    df_, remove_index = remove_outliears(df, "mw_n", 52)
    for index in remove_index:
        del signals[index]
    df_.index = wells
    # for sig, well in zip(signals, wells):
    #     sig.name = well
    #
    # df_.to_csv("computed_data.csv")
    # for sig in signals:
    #     sig.result.to_csv(f"{sig.name}.csv")

    import well_plate
    wp = well_plate.WellPlate(384, "rect")
    wp.add_data(df_["mw_n"])
    wp.plot(key="mw_n")
    print("mean :", np.mean(df_["mw_n"]))
    print("std: ", np.std(df_["mw_n"]))

    fig = px.histogram(df, x="mw_n", nbins=12)
    # fig.show()
    fig.write_html("temp.html", auto_open=True)

    print("hi")


def remove_outliears(df, column, length):
    remove_index = []
    for i in range(len(df.index)-length):
        mean = np.mean(df[column])
        errors = np.abs(df[column].to_numpy() - mean)
        max_index = np.argmax(errors)
        df = df.drop(df.index[max_index])
        remove_index.append(max_index)

    return df, remove_index


if __name__ == '__main__':
    main()

