import glob
import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)

wells = [f"{chr(64+num)}12" for num in range(1, 17)] + [f"H{num}" for num in range(1, 25)] \
        + ["A1", "A2", "B1", "B2", "A24", "P1", "P24"]


def main():
    # load data
    file_name = r"G:\Other computers\My Laptop\post_doc_2022\Data\Instrument\polymerization\RAFT8.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    for col in list(df.columns):
        if str(col).startswith("Unnamed:"):
            del df[col]

    # calibration
    def cal_func_RI(time: np.ndarray):
        return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    # load and process signals
    signals = []
    for col in df.columns:
        signals.append(chem.SECSignal(name=col, ser=df[col], cal=cal_RI,
                                      x_label="retention time (min)", y_label="signal"))

    count = 0
    df = None
    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7)
        sig.peak_picking(lb=10.6, ub=12.2)  # lb=10.85, ub=12.2

        if count == 0:
            count += 1
            df = sig.stats_df()
            sig.plot()
        else:
            df = pd.concat([df, sig.stats_df()], axis=1)

    df = df.T
    df.index = wells

    # export data
    # for sig, well in zip(signals, wells):
    #     sig.name = well
    #
    # df_.to_csv("computed_data.csv")
    # for sig in signals:
    #     sig.result.to_csv(f"{sig.name}.csv")

    import well_plate
    wp = well_plate.WellPlate(384, "rect")
    wp.add_data(df["mw_n"])
    wp.plot(key="mw_n")
    print("mean :", np.mean(df["mw_n"]))
    print("std: ", np.std(df["mw_n"]))

    # fig = px.histogram(df, x="mw_n", nbins=12)
    # fig.show()
    # fig.write_html("temp.html", auto_open=True)

    print("hi")


def main2():
    import pandas as pd
    import well_plate

    df = pd.read_csv(
        r"C:\Users\nicep\Desktop\Reseach_Post\python\chem_analysis\dev\RAFT8_data.csv", index_col=0)

    wp = well_plate.WellPlate(384, "rect")
    wp.add_data(df["mw_n"])
    wp.plot(key="mw_n")


if __name__ == '__main__':
    main2()
