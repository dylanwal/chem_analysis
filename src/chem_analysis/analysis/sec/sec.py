from functools import wraps

import numpy as np
import plotly.graph_objs as go

from src.chem_analysis.analysis.base_obj.peak import Peak
from src.chem_analysis.analysis.base_obj.chromatogram import Chromatogram


class SECPeak(Peak):
    @wraps(Peak.__init__)
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mn = None
        self.mw = None
        self.d = None
        self.std = None
        self.skew = None


class SECChromo(Chromatogram):

    def __int__(self):
        self.UV = UV
        self.RI = RI
        self.visc = visc
        self.LS

        self.cal = cal

        self.mn = None
        self.mw = None
        self.d = None
        self


    def cal_Mn_D_from_wi(MW, wi):
        # calculate Mn and D from wi vs MW data (MW goes low to high)
        data_points = len(MW)
        if MW[1] > MW[-1]:
            MW = np.flip(MW)
            wi = np.flip(wi)

        wi_d_mi = np.zeros(data_points)
        wi_m_mi = np.zeros(data_points)
        for i in range(data_points):
            if MW[i] != 0:
                wi_d_mi[i] = wi[i] / MW[i]
            wi_m_mi[i] = wi[i] * MW[i]

        Mn = np.sum(wi) / np.sum(wi_d_mi)
        Mw = np.sum(wi_m_mi) / np.sum(wi)
        D = Mw / Mn
        return Mn, D

