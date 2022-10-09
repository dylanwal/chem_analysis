import logging

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)


def main():
    # load data
    file_name = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Instrument\polymerization\RAFT4-10min-single.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)

    def cal_func_RI(time: np.ndarray):
            return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    sig = chem.SECSignal(name="sig", ser=df, cal=cal_RI, x_label="retention time (min)", y_label="RI signal")
    sig.pipeline.add(chem_bc.adaptive_polynomial_baseline)
    sig.peak_picking(lb=10.1, ub=11.8)   # lb=10.85, ub=12.2
    sig.stats()
    sig.plot(auto_open=True)

    print("hi")


if __name__ == '__main__':
    main()
