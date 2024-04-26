import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak import PeakDiscrete


class MSSignal(Signal):
    """
    Mass spectrum Signal
    """
    _peak_type = PeakDiscrete

    def __init__(self,
                 x_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "counts"
        super().__init__(x_raw, data_raw, x_label, y_label, name, id_)
