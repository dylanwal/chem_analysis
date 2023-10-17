import numpy as np

from chem_analysis.base_obj.signal_array import SignalArray
from chem_analysis.ir.ir_signal import IRSignal


class IRArray(SignalArray):
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = "wave_number",
                 y_label: str = "time",
                 z_label: str = "Absorbance",
                 name: str = None
                 ):
        super().__init__(x, y, z, x_label, y_label, z_label, name)

    # @property
    # def names(self):
    #     return [i.name for i in self.signals]
    #
    # @property
    # def y_labels(self):
    #     return [i.y_label for i in self.signals]

    @property
    def number_of_signals(self):
        return len(self.y)

    def get_signal(self, index: int) -> IRSignal:
        sig = IRSignal(x=self.x_raw, y=self.z_raw[index, :], var=self.y_raw[index], x_label=self.x_label,
                     y_label=self.z_label)
        sig.id_ = index
        sig.processor = self.processor
        return sig
