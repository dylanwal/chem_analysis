import logging
import numpy as np

from chem_analysis.base_obj.signal_discrete import SignalDiscrete
from chem_analysis.mass_spec.ms_parameters import MSParameters

logger = logging.getLogger(__name__)


class MSSignal(SignalDiscrete):
    """
    Mass spectrum Signal
    """

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
    def mz(self) -> np.ndarray:
        return self.x

    @property
    def count(self) -> np.ndarray:
        return self.y

    @property
    def high_mz(self) -> int | float:
        return np.max(np.nonzero(self.y))

    @property
    def low_mz(self) -> int | float:
        return np.min(np.nonzero(self.y))

    @property
    def number_of_peaks(self) -> int:
        return len(np.nonzero(self.y))

    @property
    def total_count(self) -> int:
        return int(np.sum(self.y))
