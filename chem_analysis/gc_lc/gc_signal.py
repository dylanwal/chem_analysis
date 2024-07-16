import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.mass_spec.ms_parameters import MSParameters


class GCSignal(Signal):
    """
    Gas Chromatogram Signal
    """
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "intensity"
        super().__init__(x, y, x_label, y_label, name, id_)
        self.parameters = parameters


