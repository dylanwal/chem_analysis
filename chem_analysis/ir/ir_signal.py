
import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.ir.ir_peak import IRPeak


class IRSignal(Signal):
    """ SECSignal
    """
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = "wave_number",
                 y_label: str = "absorbance",
                 name: str = None,
                 raw_data: bool = True,
                 _parent=None
                 ):
        super().__init__(x, y, x_label, y_label, name, raw_data, _parent)

    def __repr__(self):
        text = super().__repr__()
        return text

    @property
    def peaks(self) -> list[IRPeak]:
        return self._peaks

    def add_peak(self, peak: IRPeak):
        super().add_peak(peak)

    def delete_peak(self, peak: IRPeak | int):
        super().delete_peak(peak)

    @property
    def cm_1(self) -> np.ndarray:
        return self.x

    @property
    def micrometer(self) -> np.ndarray:
        return 1 / self.x

    @property
    def absorbance(self) -> np.ndarray:
        return self.y

    @property
    def transmittance(self) -> np.ndarray:
        return np.exp(-self.y)
