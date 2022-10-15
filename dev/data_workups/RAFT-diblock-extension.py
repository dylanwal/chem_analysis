import glob
import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)

def main():
    # load data
    file_name = r"G:\Other computers\My Laptop\post_doc_2022\Data\Instrument\polymerization\homo_diblock_data.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)

    # calibration
    def cal_func_RI(time: np.ndarray):
        return 10 ** (-0.6 * time + 10.644)

    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=1_090_000, name="RI calibration")

    block1 = chem.SECSignal(name="block 1", ser=df[df.columns[0]], cal=cal_RI,
                           x_label="retention time (min)", y_label="signal")

    diblock = chem.SECSignal(name="diblock", ser=df[df.columns[1]], cal=cal_RI,
                           x_label="retention time (min)", y_label="signal")

    block1.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7)
    block1.peak_picking(lb=10, ub=12)
    print(block1.stats_df())
    fig = block1.plot(auto_open=False, op_cal=False, normalize=True)

    diblock.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7)
    diblock.peak_picking(lb=9.4, ub=11.5)
    print(diblock.stats_df())
    fig = diblock.plot(fig, normalize=True, auto_open=False)
    fig.update_layout(legend=dict(x=.5, y=.95))
    fig.data[4].fillcolor = 'rgba(172,24,25,0.5)'
    fig.write_html('temp.html', auto_open=True)

    print("hi")


if __name__ == '__main__':
    main()
