from enum import Enum

import numpy as np

from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak_SEC import PeakSEC


class SECTypes(Enum):
    UNKNOWN = -1
    RI = 0
    UV = 1
    LS = 2
    VISC = 3


class SECSignal(Signal):
    """
    SECSignal
    """
    TYPES_ = SECTypes
    _peak_type = PeakSEC

    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 calibration: SECCalibration = None,
                 type_: SECTypes = SECTypes.UNKNOWN,
                 x_label: str = "retention time",
                 y_label: str = "signal",
                 name: str = None,
                 id_: int = None
                 ):
        super().__init__(x_raw, y_raw, x_label, y_label, name, id_)
        self.calibration = calibration
        self.type_ = type_

    def __repr__(self):
        text = super().__repr__()
        if self.type_ is not SECTypes.UNKNOWN:
            text += f" ({self.type_.name})"

        return text

    @property
    def mw_i(self) -> np.ndarray | None:
        if self.calibration is None:
            return None

        return self.calibration.get_y(self.x)
