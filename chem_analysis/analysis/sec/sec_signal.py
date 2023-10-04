from enum import Enum

import numpy as np

from chem_analysis.analysis.base_obj.calibration import Calibration
from chem_analysis.analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.sec.sec_peak import SECPeak


class SECType(Enum):
    RI = 0
    UV = 1
    LS = 2
    VISC = 3


class SECSignal(Signal):
    """ SECSignal

    Extends Signal for SEC (size exclusion chronograph)

    Attributes
    ----------
    cal: Cal
        calibration
    peaks: List[SECPeak]


    """
    _peak = SECPeak

    def __init__(self, x: np.ndarray,
                 y: np.ndarray,
                 calibration: Calibration,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 raw_data: bool = True,
                 _parent=None
                 ):
        super().__init__(x, y, x_label, y_label, name, raw_data, _parent)
        self.calibration = calibration

    def mw_i(self) -> np.ndarray:
        return self.calibration(self.result.index.to_numpy())

    def calc(self):
        super().calc()

        mol_weight = self.cal(self.result.index.to_numpy())
        self._result_weight = pd.Series(self.result.to_numpy(), mol_weight)
