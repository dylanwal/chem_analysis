import numpy as np
import pywt

from chem_analysis.processing.processing_method import Baseline


def wavelet_baseline_correction(y, wavelet='sym6', level=6):
    coeffs = pywt.wavedec(y, wavelet, mode='smooth', level=level)

    # Zero out detail coefficients to retain only the approximation (baseline)
    coeffs[1:] = [np.zeros_like(c) for c in coeffs[1:]]
    baseline = pywt.waverec(coeffs, wavelet, mode='smooth')
    baseline = baseline[:len(y)]
    return baseline


class Wavelet(Baseline):
    def __init__(self,
                 wavelet: str = 'sym6',
                 level: int = 6,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.wavelet = wavelet
        self.level = level

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return wavelet_baseline_correction(y, self.wavelet, self.level)
