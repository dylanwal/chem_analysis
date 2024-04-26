
import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.analysis.peak import PeakDiscrete


class IRSignal2D(Signal2D):
    _signal = MSSignal
    _peak_type = PeakDiscrete

    def __init__(self,
                 x_raw: np.ndarray,
                 time_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "time"
        z_label = z_label or "counts"
        super().__init__(x_raw, time_raw, data_raw, x_label, y_label, z_label, name)

    def get_signal(self, index: int, processed: bool = False) -> MSSignal:
        return super().get_signal(index, processed)
