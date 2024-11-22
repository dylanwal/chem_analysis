from __future__ import annotations
import dataclasses
import logging

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.sec.sec_math_functions import calculate_Mn_D_from_wi
from chem_analysis.analysis.peak import PeakParent, PeakContinuous
from chem_analysis.sec.sec_calibration import SECCalibration

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class PeakParentSEC(PeakParent):
    mw_i: np.ndarray | None
    calibration: SECCalibration | None


class PeakSEC(PeakContinuous):
    def __init__(self, parent: PeakParentSEC, x: np.ndarray, y: np.ndarray, id_: int = None):
        super().__init__(parent, x, y, id_)
        self.parent: PeakParentSEC = parent  # duplicate, but helps with type hinting

        self._mw_n = None
        self._mw_d = None
        self._w_i = None
        self._x_i = None
        self._mw_i_None = False

    def _get_exclude_from_stats(self) -> set[str]:
        return super()._get_exclude_from_stats().union({"x_i", "w_i", "mw_i"})

    @property
    def mw_i(self) -> np.ndarray | None:
        """ molecular weight of i-mer """
        if self.parent.mw_i is None:
            return None
        if self._mw_i_None:
            return None
        mw_i = self.parent.mw_i[self.bounds]
        if np.min(mw_i) == 0 and np.max(mw_i) == 0:
            # if mw_i is all zero, then peak is outside calibration and we shouldn't try to do calculations around it
            self._mw_i_None = True
            # self.stats = PeakStats(self)
            return None

        return mw_i

    @property
    def w_i(self) -> np.ndarray | None:
        """ w_i = mole fraction , weight fraction of i-mer"""
        if self.mw_i is None:
            return None
        if self._w_i is None:
            self._w_i = self.y / np.abs(np.trapz(x=self.mw_i, y=self.y))  # abs because mw_i and y are reversed order
        return self._w_i

    @property
    def mw_n(self) -> float | None:
        """ number average molecular weight """
        if self.mw_i is None:
            return None
        if self._mw_n is None:
            self._mw_n, self._mw_d = calculate_Mn_D_from_wi(mw_i=self.mw_i, wi=self.w_i)
        return self._mw_n

    @property
    def mw_d(self) -> float | None:
        """ dispersity of molecular weight """
        if self.mw_i is None:
            return None
        if self._mw_n is None:
            self._mw_n, self._mw_d = calculate_Mn_D_from_wi(mw_i=self.mw_i, wi=self.w_i)
        return self._mw_d

    @property
    def mw_w(self) -> float | None:
        """ weight average molecular weight """
        if self.mw_i is None:
            return None
        return self.mw_n * self.mw_d

    @property
    def x_i(self) -> np.ndarray | None:
        """ x_i = mole fraction,   mole fraction of i-mer """
        if self.mw_i is None:
            return None
        if self._x_i is None:
            numerator = self.w_i * self.mw_n
            zero_mask = (self.mw_i == 0)
            if np.sum(zero_mask) != 0:
                logger.warning("PeakContinuous extends outside calibration window. Results may contain errors.")
            self._x_i = np.divide(numerator, self.mw_i, out=np.zeros_like(numerator), where=~zero_mask)
            # self._x_i = self.w_i * self.mw_n / self.mw_i

        return self._x_i

    @property
    def mw_max(self) -> float:
        return np.max(self.mw_i)

    @property
    def mw_mean(self) -> float:
        return general_math.get_mean_of_pdf(self.mw_i, y_norm=self.x_i)

    @property
    def mw_std(self):
        return general_math.get_standard_deviation_of_pdf(self.mw_i, y_norm=self.x_i, mean=self.mw_mean)

    @property
    def mw_skew(self):
        return general_math.get_skew_of_pdf(self.mw_i, y_norm=self.x_i, mean=self.mw_mean,
                                            standard_deviation=self.mw_std)

    @property
    def mw_kurtosis(self):
        return general_math.get_kurtosis_of_pdf(self.mw_i, y_norm=self.x_i, mean=self.mw_mean,
                                                standard_deviation=self.mw_std)

    @property
    def mw_fwhm(self):
        """mw_full_width_half_max"""
        return general_math.get_full_width_at_height(x=self.mw_i, y=self.x_i, height=0.5)

    @property
    def mw_asym(self):
        """mw_asymmetry_factor"""
        return general_math.get_asymmetry_factor(x=self.mw_i, y=self.x_i, height=0.1)
