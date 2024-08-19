from __future__ import annotations

import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.sec.sec_calibration import SECCalibration
from chem_analysis.sec.sec_signal import SECSignal, SECTypes
from chem_analysis.analysis.integration.sec_peak import PeakSEC


class SECSignalArray(Signal2D):
    TYPES_ = SECTypes
    _signal = SECSignal
    _peak_type = PeakSEC

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 calibration: SECCalibration = None,
                 type_: SECTypes = SECTypes.UNKNOWN,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "retention_time"
        y_label = y_label or "time"
        z_label = z_label or "signal"
        super().__init__(x, y, z, x_label, y_label, z_label, name)
        self.calibration = calibration
        self.type_ = type_

    def get_signal(self, index: int, processed: bool = True) -> SECSignal:
        sig = super().get_signal(index, processed)
        sig.calibration = self.calibration
        return sig

    # @classmethod
    # def from_file(cls, path: str | pathlib.Path, calibration: SECCalibration = None) -> SECSignalArray:
    #     class_ = super().from_file(path)
    #     class_.calibration = calibration
    #
    #     return class_
