
import numpy as np

from chem_analysis.base_obj.signal_discrete_2d import SignalDiscrete2D
from chem_analysis.mass_spec.ms_signal_low_res import MSSignalLowRes


class MSSignalLowRes2D(SignalDiscrete2D):
    _signal = MSSignalLowRes

    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "time"
        z_label = z_label or "counts"
        super().__init__(x_raw, y_raw, data_raw, x_label, y_label, z_label, name, id_)

    def get_signal(self, index: int, processed: bool = False) -> MSSignalLowRes:
        return super().get_signal(index, processed)
