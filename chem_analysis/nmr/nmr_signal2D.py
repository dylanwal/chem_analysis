import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.nmr.nmr_signal import NMRSignal
from chem_analysis.nmr.NMR_parameters import NMRParameters


class NMRSignal2D(Signal2D):
    _signal = NMRSignal

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 parameters: NMRParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "ppm"
        y_label = y_label or "time"
        z_label = z_label or "signal"
        super().__init__(x, y, z, x_label, y_label, z_label, name, id_)
        self.parameters = parameters

    def get_signal(self, index: int, processed: bool = True) -> NMRSignal:
        return super().get_signal(index, processed)
