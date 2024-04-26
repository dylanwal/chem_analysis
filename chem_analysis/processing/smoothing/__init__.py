from chem_analysis.processing.smoothing.savitzky_golay import SavitzkyGolay
from chem_analysis.processing.smoothing.gaussian import Gaussian
from chem_analysis.processing.smoothing.time_smoothing import ExponentialTime, GaussianTime
from chem_analysis.processing.smoothing.rolling_window import RollingWindow

# TODO:
# Denoise
# non-local means: noise factor 0.75 blockwise
# median-modified wiener
# gaussian


# Savitzky-Golay
# moving average: span
# Whittaker Smoother: smooth factor 36
# Wavelets: scales=4, fraction=1%
# Cadzow


# class LineBroadening(Smoothing):
#     def __init__(self, degree: float = 1):
#         self.degree = degree  # Hz
#         self._y_baseline = None
#
#     @property
#     def y_baseline(self) -> np.ndarray:
#         return self._y_baseline
#
#     def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
#         length = len(nmrData.allFid[-1][0])
#         sp.multiply(nmrData.allFid[-1][:], sp.exp(-nmrData.fidTimeForLB[:length] * self.degree * np.pi))
#         return x, y
