import numpy as np

from chem_analysis.analysis.base_obj.calibration import Calibration
from chem_analysis.analysis.base_obj.peak import Peak, PeakSupports


class SECPeakSupports(PeakSupports):
    cal: Calibration


class SECPeak(Peak):
    """ SECPeak

    Extends the Peak class to for SEC (size exclusion chromatograph).

    Attributes
    ----------
    mw_i: np.ndarray
        molecular weight of i-mer
    wi: np.ndarray
        weight fraction of i-mer
    xi: np.ndarray
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
        self.mw_i = None
        self.wi = None
        self.xi = None

        self.mw_n = None
        self.mw_w = None
        self.mw_max = None
        self.mw_d = None
        self.mw_mean = None
        self.mw_std = None
        self.mw_skew = None
        self.mw_kurtosis = None
        self.mw_full_width_half_max = None
        self.mw_asymmetry_factor = None

        super().__init__(parent, low_bound_index, high_bound_index, id_)

    def __repr__(self):
        return super().__repr__() + f" | Mn: {self.mw_n} (D: {self.mw_d})"

    def calculate_parameters(self):
        """ Calculates the molecular weight stats for the peak. """
        super().calculate_parameters()

        mw_i = np.flip(self._parent.result_weight.index[self.low].to_numpy())
        # 'np.flip' so small mw is first avoids having to flip everything in calc below
        self.wi = np.flip(self._parent.result_weight.iloc[self.slice].to_numpy() / \
                  np.trapz(x=self.mw_i, y=self._parent.result_weight.iloc[self.slice].to_numpy()))

        self.mw_n, self.mw_d = cal_Mn_D_from_wi(mw_i=self.mw_i, wi=self.wi)
        self.mw_w = self.mw_n * self.mw_d
        self.mw_max = self._parent.cal(self.max_loc)

        self.xi = self.wi * self.mw_n / self.mw_i

        self.mw_mean = np.trapz(x=self.mw_i, y=self.mw_i*self.xi)
        self.mw_std = np.sqrt(np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 2))
        self.mw_skew = np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 3) / self.mw_std ** 3
        self.mw_kurtosis = (np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 4) / self.mw_std ** 4)- 3
        self.mw_full_width_half_max = self.get_fw(x=self.mw_i, y=self.xi, height=0.5)
        self.mw_asymmetry_factor = self.get_asym(x=self.mw_i, y=self.xi, height=0.1)

    def stats(self, op_print: bool = True, op_headers: bool = True, window: int = 150, headers: dict = None):
        if headers is None:
            headers = {  # attribute: print
                "id_": "id", "lb_loc": "low bound", "max_loc": "max", "hb_loc": "high bound", "area": "area",
                "mw_n": "mw_n", "mw_w": "mw_w", "mw_d": "mw_d", "mw_max": "mw_max", "mw_std": "mw_std",
                "mw_skew": "mw_skew", "mw_kurtosis": "mw_kurtosis", "mw_asym": "mw_asym"
            }

        return super().stats(op_print, op_headers, window, headers)
