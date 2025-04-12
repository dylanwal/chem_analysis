import logging
import pathlib

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

    def update_mz(self, mz: np.ndarray):
        if len(mz.shape) != 1:
            raise ValueError("Mass spectrum must have exactly one mass spectrum")
        if mz[0] > mz[-1]:
            mz = np.flip(mz)

        indexes = np.nonzero(self.y != 0)[0]
        if np.max(mz) < self.mz[indexes[-1]] or np.min(mz) > self.mz[indexes[0]]:
            logger.warning("The new 'mz' is smaller than the old 'mz' and some signal may be lost.")

        y = np.zeros_like(mz, dtype=self.y.dtype)
        min_, max_ = np.min(mz), np.max(mz)
        for i, mz_ in enumerate(self.mz):
            if min_ <= mz_ <= max_:
                index = np.nonzero(mz == mz_)[0][0]
                y[index] = self.y[i]

        self.y = y
        self.x = mz

    def write_csv(self,
               path: str | pathlib.Path,
               headers: bool = False,
               encoding: str = "utf-8",
               **kwargs
               ):
        if "fmt" not in kwargs:
            kwargs["fmt"] = ("%.0i", "%.0f")

        super().write_csv(path, headers, encoding, **kwargs)
