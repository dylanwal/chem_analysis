import logging
import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.mass_spec.ms_parameters import MSParameters

logger = logging.getLogger(__name__)


class MSSignal(Signal):
    """
    Mass spectrum Signal
    """
    _discrete = True

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "counts"
        super().__init__(x, y, x_label, y_label, name, id_)
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

    def get_intensity(self, mz: int | float) -> int | float:
        if mz in self.mz:
            index = np.where(self.mz == mz)[0][0]
            return self.y[index]

        return 0
