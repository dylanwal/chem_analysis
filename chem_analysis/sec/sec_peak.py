from collections import OrderedDict

import numpy as np

from chem_analysis.base_obj.calibration import Calibration
from chem_analysis.base_obj.peak import Peak, PeakSupports
from chem_analysis.sec.sec_math_functions import calculate_Mn_D_from_wi
import chem_analysis.algorithms.general_math as general_math


class SECPeakSupports(PeakSupports):
    calibration: Calibration
    mw_i: np.ndarray | None


class SECPeak(Peak):
    """ SECPeak

    Extends the Peak class to for SEC (size exclusion chromatograph).

    Attributes
    ----------
    mw_i: np.ndarray
        molecular weight of i-mer
    w_i: np.ndarray
        weight fraction of i-mer
    x_i: np.ndarray
        mole fraction of i-mer
    mw_n: float
        number average molecular weight
    mw_w: float
        weight average molecular weight
    mw_max: float
        molecular weight at peak max
    mw_d: float
        dispersity of molecular weight
    mw_mean: float
        mean molecular weight distribution (same as mw_n)
    mw_std: float
        standard deviation of molecular weight distribution
    mw_skew: float
        skew of molecular weight distribution
        symmetric: -0.5 to 0.5; moderate skew: -1 to -0.5 or 0.5 to 1; high skew: <-1 or >1;
        positive tailing to higher numbers; negative tailing to smaller numbers
    mw_kurtosis: float
        kurtosis of molecular weight distribution (Fisher)  (Warning: highly sensitive to peak bounds)
        negative: flatter peak; positive: sharp peak
    mw_full_width_half_max: float
        full width half max of molecular weight distribution
    mw_asymmetry_factor: float
        asymmetry factor of molecular weight distribution; distance from the center line of the peak to the back
        slope divided by the distance from the center line of the peak to the front slope
        >1 tailing to larger values; <1 tailing to smaller numbers

    """
    __slots__ = ("mw_max", "mw_mean", "mw_std", "mw_skew", "mw_kurtosis", "mw_full_width_half_max",
                 "mw_asymmetry_factor")
    _stats = ["mw_n", "mw_w", "mw_d", "mw_max", "mw_std", "mw_skew", "mw_kurtosis", "mw_full_width_half_max",
              "mw_asymmetry_factor"]

    def __init__(self, parent: SECPeakSupports, low_bound_index: int, high_bound_index: int, id_: int = None):
        """
        Parameters
        ----------
        parent: SECPeakSupports
            parent object that the peak is associated with
        low_bound_index: int
            lower bound index
        high_bound_index: int
            higher bound index

        """
        self.mw_max = None
        self.mw_mean = None
        self.mw_std = None
        self.mw_skew = None
        self.mw_kurtosis = None
        self.mw_full_width_half_max = None
        self.mw_asymmetry_factor = None

        super().__init__(parent, low_bound_index, high_bound_index, id_)

    def __repr__(self):
        return super().__repr__() + f" | Mn: {self.mw_n} (D: {self.mw_d})"

    @property
    def mw_i(self) -> np.ndarray | None:
        if self.parent.mw_i is None:
            return None
        return self.parent.mw_i[self.low_bound_index:self.high_bound_index]

    @property
    def w_i(self) -> np.ndarray | None:
        """ w_i = mole fraction """
        if self.mw_i is None:
            return None
        return self.y / np.trapz(x=self.mw_i, y=self.y)

    @property
    def mw_n(self) -> float | None:
        if self.mw_i is None:
            return None
        mw_n, mw_d = calculate_Mn_D_from_wi(mw_i=self.mw_i, wi=self.w_i)
        return mw_n

    @property
    def mw_d(self) -> float | None:
        if self.mw_i is None:
            return None
        mw_n, mw_d = calculate_Mn_D_from_wi(mw_i=self.mw_i, wi=self.w_i)
        return mw_d

    @property
    def mw_w(self) -> float | None:
        if self.mw_i is None:
            return None
        return self.mw_n * self.mw_d

    @property
    def x_i(self) -> np.ndarray | None:
        """ x_i = mole fraction """
        if self.mw_i is None:
            return None

        return self.w_i * self.mw_n / self.mw_i

    def calculate_parameters(self):
        """ Calculates the molecular weight stats for the peak. """
        super().calculate_parameters()

        if self.mw_i is None:
            # no calibration present
            return

        mw_i = np.flip(self.mw_i)
        x_i = np.flip(self.x_i)
        # 'np.flip' so small mw is first avoids having to flip everything in calc below

        self.mw_max = np.argmax(x_i) + self.low_bound_index
        self.mw_mean = general_math.get_mean_of_pdf(mw_i, y_norm=x_i)
        self.mw_std = general_math.get_standard_deviation_of_pdf(mw_i, y_norm=x_i, mean=self.mw_mean)
        self.mw_skew = general_math.get_skew_of_pdf(mw_i, y_norm=x_i, mean=self.mw_mean,
                                                    standard_deviation=self.mw_std)
        self.mw_kurtosis = general_math.get_kurtosis_of_pdf(mw_i, y_norm=x_i, mean=self.mw_mean,
                                                            standard_deviation=self.mw_std)
        self.mw_full_width_half_max = general_math.get_full_width_at_height(x=mw_i, y=x_i, height=0.5)
        self.mw_asymmetry_factor = general_math.get_asymmetry_factor(x=mw_i, y=x_i, height=0.1)

    def stats(self) -> OrderedDict:
        dict_ = super().stats()
        for stat in self._stats:
            dict_[stat] = getattr(self, stat)
        return dict_
