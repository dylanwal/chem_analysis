
import numpy as np

from chem_analysis.base_obj.signal_array import SignalArray
from chem_analysis.nmr.nmr_signal import NMRSignal
from chem_analysis.analysis.peak import PeakBoundedStats


class NMRSignalArray(SignalArray):
    _signal = NMRSignal
    _peak_type = PeakBoundedStats

    def __init__(self,
                 x_raw: np.ndarray,
                 time_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "ppm"
        y_label = y_label or "time"
        z_label = z_label or "signal"
        super().__init__(x_raw, time_raw, data_raw, x_label, y_label, z_label, name)

    def get_signal(self, index: int) -> NMRSignal:
        sig = NMRSignal(x_raw=self.x_raw, y_raw=self.raw_data[index, :], x_label=self.x_label, y_label=self.y_label,
                       name=f"time: {self.time[index]}", id_=index)
        sig.time = self.time[index]
        sig.processor = self.processor.get_copy()
        return sig
