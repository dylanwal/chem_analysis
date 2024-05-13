from __future__ import annotations
import numpy as np

from chem_analysis.base_obj.signal_discrete import SignalDiscrete
from chem_analysis.mass_spec.ms_parameters import MSParameters


class MSSignalLowRes(SignalDiscrete):
    """
    Mass spectrum Signal
    Low resolution (LowRes) limits mz values to integers.
    """

    def __init__(self,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "counts"
        super().__init__(np.arange(len(data_raw)), data_raw, x_label, y_label, name, id_)
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

    @classmethod
    def from_peaks(cls, mz: np.ndarray, data_raw: np.ndarray) -> MSSignalLowRes:
        data_new = np.zeros(np.max(mz))
        mz = np.round(mz).astype('uint64')
        data_new[mz] = data_raw
        return cls(data_new)
