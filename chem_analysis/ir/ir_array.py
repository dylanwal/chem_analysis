
import numpy as np

from chem_analysis.base_obj.signal_array import SignalArray
from chem_analysis.ir.ir_signal import IRSignal
from chem_analysis.analysis.peak_IR import PeakIR


class IRSignalArray(SignalArray):
    _peak_type = PeakIR

    def __init__(self,
                 x_raw: np.ndarray,
                 time_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = "wave_number",
                 y_label: str = "time",
                 z_label: str = "absorbance",
                 name: str = None
                 ):
        super().__init__(x_raw, time_raw, data_raw, x_label, y_label, z_label, name)

    @property
    def cm_1(self) -> np.ndarray:
        return self.x

    @property
    def micrometer(self) -> np.ndarray:
        return 1 / self.x * 1000

    @property
    def absorbance(self) -> np.ndarray:
        return self.data

    @property
    def transmittance(self) -> np.ndarray:
        return np.exp(-self.data)

    def get_signal(self, index: int) -> IRSignal:
        sig = IRSignal(x_raw=self.x, y_raw=self.data[index, :], x_label=self.x_label, y_label=self.y_label,
                       name=f"time: {self.time[index]}", id_=index)
        sig.processor = self.processor
        return sig
