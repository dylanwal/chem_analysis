import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)


def extract_from_excel2(file_name: str):
    return pd.read_excel(file_name, sheet_name=None, names=["time_RI", "RI", "timeUV", "UV"], skiprows=[1])


def main():
    file_name = r"G:\Other computers\My Laptop\post_doc_2022\Data\Instrument\polymerization\DW-RAFT-2-SEC.xlsx"
    data = extract_from_excel2(file_name)

    def cal_func_RI(time: np.ndarray):
        return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    signals = []
    for name, df in data.items():
        signals.append(chem.SECSignal(name=name, xy=df[["time_RI", "RI"]].to_numpy(), cal=cal_RI))

    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline)
        sig.peak_picking(lb=10, ub=11.8)
        sig.stats(num_sig_figs=4, op_headers=False)
        sig.plot()

    print("hi")


if __name__ == '__main__':
    main()
