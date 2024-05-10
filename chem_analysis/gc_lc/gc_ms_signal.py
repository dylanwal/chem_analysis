import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.gc_lc.gc_parameters import GCParameters


class GCMSSignal(Signal):
    """
    Gas Chromatogram Signal
    """
    def __init__(self,
                 x_raw: np.ndarray,
                 ms_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: GCParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "intensity"
        self.ms_raw = ms_raw
        data_raw = np.sum(ms_raw, axis=1)
        super().__init__(x_raw, data_raw, x_label, y_label, name, id_)
        self.parameters = parameters

    def get_ms(self, index: int | slice):
        ...
