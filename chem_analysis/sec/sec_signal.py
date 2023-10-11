from enum import Enum

import numpy as np

from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.sec.sec_peak import SECPeak


class SECTypes(Enum):
    UNKNOWN = -1
    RI = 0
    UV = 1
    LS = 2
    VISC = 3


class SECSignal(Signal):
    """ SECSignal
    """
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 calibration: SECCalibration = None,
                 type_: SECTypes = SECTypes.UNKNOWN,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 raw_data: bool = True,
                 _parent=None
                 ):
        super().__init__(x, y, x_label, y_label, name, raw_data, _parent)
        self.calibration = calibration
        self.type_ = type_

    def __repr__(self):
        text = super().__repr__()
        if self.type_ is not SECTypes.UNKNOWN:
            text += f" ({self.type_.name})"

        return text

    @property
    def peaks(self) -> list[SECPeak]:
        return self._peaks

    def add_peak(self, peak: SECPeak):
        super().add_peak(peak)

    def delete_peak(self, peak: SECPeak | int):
        super().delete_peak(peak)

    @property
    def mw_i(self) -> np.ndarray | None:
        if self.calibration is None:
            return None

        return self.calibration.get_y(self.x)
