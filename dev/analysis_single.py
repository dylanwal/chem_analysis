import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)


def main():
    # load data
    file_name = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Polyester\publish\SEC\DW1-3-1-redo[DW1-8].csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    df = df.iloc[:1772, :10]
    df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    df.index.names = ["time (min)"]
    df = df.apply(pd.to_numeric, errors='coerce')

    # define calibration
    def cal_func_RI(time: np.ndarray):
        return 10 ** (0.0206 * time ** 2 - 1.0422 * time + 15.081)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    # define signal
    signal = chem.SECSignal(ser=df["RI"], cal=cal_RI, x_label="retention time (min)", y_label="signal")

    # data processing
    signal.pipeline.add(chem_bc.adaptive_polynomial_baseline)
    signal.peak_picking(lb=13, ub=18)

    # print stats / plotting
    signal.stats(num_sig_figs=4)
    signal.plot()

    print("hi")


if __name__ == '__main__':
    main()

