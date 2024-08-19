
import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.ir.ir_signal import IRSignal


class IRSignal2D(Signal2D):
    _signal = IRSignal

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "wave_number"
        y_label = y_label or "time"
        z_label = z_label or "absorbance"
        super().__init__(x, y, z, x_label, y_label, z_label, name, id_)

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

    def get_signal(self, index: int, processed: bool = True) -> IRSignal:
        return super().get_signal(index, processed)
