from typing import Union, Callable

import pandas as pd

from src.chem_analysis.analysis.base_obj.calibration import Cal
from src.chem_analysis.analysis.base_obj.chromatogram import Chromatogram
from src.chem_analysis.analysis.base_obj.signal_ import Signal
from src.chem_analysis.analysis.sec.sec_signal import SECSignal


class SECChromo(Chromatogram):

    __signal = SECSignal

    def __init__(self, data: Union[pd.DataFrame, Signal, list[Signal]], cal: Union[Cal, Callable] = None):
        super().__init__(data)
        if not isinstance(cal, Cal):
            cal = Cal(cal)
        self.cal = cal


def local_run():
    from scipy.stats import norm
    import numpy as np

    def cal(time: np.ndarray):
        return 0.0167*time**2 - 0.9225*time + 14.087

    nx = 1000
    ny = 3
    x = np.linspace(0, 25, nx)
    y = np.empty((ny, nx))
    for i in range(ny):
        rv = norm(loc=15, scale=0.6)
        rv2 = norm(loc=18, scale=0.6)
        y[i, :] = 5 * np.linspace(0, 1, nx) + np.random.random(nx) + 100 * rv.pdf(x) + 20 * rv2.pdf(x)
    df = pd.DataFrame(data=y.T, index=x)
    df.columns = ["RI", "UV", "LS"]
    df.index.names = ["time"]

    chro = SECChromo(data=df, cal=cal)
    chro.baseline(deg=1)
    chro.auto_peak_picking()
    chro.plot()
    chro.stats()
    print("done")


if __name__ == "__main__":
    local_run()

    # path = r"C:\Users\nicep\Desktop\Reseach_Post\Data\Polyester\DW1_3\SEC\DW1-3-2[DW1-3].csv"
    # df = pd.read_csv(path, header=0, index_col=0)
    # df = df.iloc[:, :10]
    # df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    # df.index.names = ["time (min)"]


