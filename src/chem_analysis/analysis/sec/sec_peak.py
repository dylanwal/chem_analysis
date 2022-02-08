from typing import Protocol, Union

import numpy as np
import pandas as pd

from src.chem_analysis.analysis.base_obj.calibration import Cal
from src.chem_analysis.analysis.base_obj.peak import Peak


def cal_Mn_D_from_wi(mw_i: np.ndarray, wi: np.ndarray) -> tuple[float, float]:
    # calculate Mn and D from wi vs MW data (MW goes low to high)
    data_points = len(mw_i)
    if mw_i[1] > mw_i[-1]:
        mw_i = np.flip(mw_i)
        wi = np.flip(wi)

    wi_d_mi = np.zeros(data_points)
    wi_m_mi = np.zeros(data_points)
    for i in range(data_points):
        if mw_i[i] != 0:
            wi_d_mi[i] = wi[i] / mw_i[i]
        wi_m_mi[i] = wi[i] * mw_i[i]

    mw_n = np.sum(wi) / np.sum(wi_d_mi)
    mw_w = np.sum(wi_m_mi) / np.sum(wi)
    mw_d = mw_w / mw_n
    return mw_n, mw_d


class SECPeakSupports(Protocol):
    name: str
    result: pd.Series
    result_weight: pd.Series
    cal: Cal


class SECPeak(Peak):

    def __init__(self, parent: SECPeakSupports, lb_index: float, hb_index: float):
        super().__init__(parent, lb_index, hb_index)

        self.mw_n = None
        self.mw_w = None
        self.mw_max = None
        self.mw_d = None
        self.mw_mean = None
        self.mw_std = None
        self.mw_skew = None
        self.mw_kurtosis = None

        self.calc()

    def __repr__(self):
        return f"peak: {self.id_} at {self.max_loc} with Mn: {self.mw_n} (D: {self.mw_d})"

    @property
    def x_weight(self) -> np.ndarray:
        return self._parent.result_weight.index[self.slice_].to_numpy()

    @property
    def y_weight(self) -> np.ndarray:
        return self._parent.result_weight.iloc[self.slice_].to_numpy()

    def calc(self):
        super().calc()

        self.mw_n, self.mw_d = cal_Mn_D_from_wi(mw_i=self.x_weight, wi=self.y_weight)
        self.mw_w = self.mw_n * self.mw_d
        self.mw_max = self._parent.cal(self.max_loc)

        wi = self.y_weight
        mw_i = self.x_weight
        xi = wi*self.mw_n/self.x_weight

        self.mw_mean = self.mw_n
        self.mw_std = np.sqrt(np.sum(xi*(mw_i-self.mw_mean)**2))
        self.mw_skew =

    def stats(self, op_print: bool = True, op_headers: bool = True, window: int = 150, headers: dict = None):
        if headers is None:
            headers = {  # attribute: print
                "id_": "id", "lb_loc": "low bound", "max_loc": "max", "hb_loc": "high bound", "area": "area",
                "mw_n": "mw_n", "mw_w": "mw_w", "mw_d": "mw_d", "mw_max": "mw_max"
            }

        return super().stats(op_print, op_headers, window, headers)
