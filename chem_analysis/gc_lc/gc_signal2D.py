import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.gc_lc.gc_parameters import GCParameters


class GCSignal2D(Signal2D):
    """
    Gas Chromatogram Signal
    """
    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 z_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 parameters: GCParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "time"
        z_label = z_label or "intensity"
        super().__init__(x_raw, y_raw, z_raw, x_label, y_label, z_label, name, id_)
        self.parameters = parameters

