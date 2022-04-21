import logging

import numpy as np
import pandas as pd

import chem_analysis.analysis as chem

chem.logger_analysis.setLevel(logging.CRITICAL)


def single_analysis(path: str):
    df = pd.read_csv(path, header=0, index_col=0)
    df = df.iloc[:1772, :10]
    df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    df.index.names = ["time (min)"]
    df = df.apply(pd.to_numeric, errors='coerce')


    def cal_func_UV(time: np.ndarray):
        return 10 ** (0.02 * time ** 2 - 1.0175 * time + 14.841)


    def cal_func_RI(time: np.ndarray):
        return 10 ** (0.0206 * time ** 2 - 1.0422 * time + 15.081)

    cal_UV = chem.Cal(cal_func_UV, lb=160, ub=319_000, name="UV calibration")
    cal_RI = chem.Cal(cal_func_RI, lb=160, ub=319_000, name="RI calibration")

    sig_RI = chem.SECSignal(ser=df["RI"], cal=cal_RI)
    sig_RI.auto_full()
    sig_RI.plot(title=path)
    sig_RI.stats()

    sig_UV = chem.SECSignal(ser=df["UV"], cal=cal_UV)
    sig_RI.auto_full()
    sig_UV.plot(title=path)
    sig_UV.stats()


def run_muti():
    import glob
    import os

    folder = r"C:\Users\nicep\Desktop\Reseach_Post\Data\Polyester\publish\SEC"
    os.chdir(folder)
    files = glob.glob("*.csv")
    for file in files[0:5]:
        print(file)
        single_analysis(file)


def run():
    # path = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Polyester\publish\SEC\DW1-3-1-redo[DW1-8].csv"
    path = r"G:\Other computers\My Laptop\post_doc_2022\Data\Polyester\publish\SEC\DW1-3-1-redo[DW1-8].csv"
    single_analysis(path)


if __name__ == '__main__':
    run()
    # run_muti()
