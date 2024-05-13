import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.mass_spec.ms_parameters import MSParameters
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
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "counts"
        super().__init__(x_raw, data_raw, x_label, y_label, name, id_)
        self.parameters = parameters

    @property
    def high_mz(self) -> int:
        return np.max(np.nonzero(self.y))

    @property
    def low_mz(self) -> int:
        return np.min(np.nonzero(self.y))

    @property
    def number_of_peaks(self) -> int:
        return len(np.nonzero(self.y))

    @property
    def total_count(self) -> int:
        return int(np.sum(self.y))
