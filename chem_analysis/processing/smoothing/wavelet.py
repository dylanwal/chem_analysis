
import numpy as np


from chem_analysis.processing.smoothing.savitzky_golay import Smoothing


class Wavelet(Smoothing):
    """
    Steps for Wavelet Denoising
        1) Decompose the Signal:
            Perform a Discrete Wavelet Transform (DWT) to break the signal into approximation (low-frequency) and detail (high-frequency) components at multiple levels.

        2) Threshold the Details:
            Apply a thresholding technique to the detail coefficients to suppress noise while preserving important signal features.
            Two common thresholding strategies:
            Hard Thresholding: Set all coefficients below a threshold to zero.
            Soft Thresholding: Shrink coefficients toward zero, making the transition smoother.

        3) Reconstruct the Signal:
            Perform the inverse DWT using the thresholded coefficients to reconstruct the denoised signal.

    Why Use Wavelets for Denoising?
    * Multi-resolution: Wavelets can separate the noise at different scales from the true signal, which is difficult with other techniques like Fourier transforms.
    * Localized Analysis: Wavelets work well for signals with non-stationary features (e.g., sudden spikes, edges).
    * Efficient: Compared to traditional filtering techniques, wavelet-based denoising can achieve better results with less computational overhead.
    """
    def __init__(self, wavelet: str = 'db4', threshold: float = 0.5, level: int = 4, temporal_processing: int = 1):
        """
        https://pywavelets.readthedocs.io/en/latest/regression/wavelet.html

        Parameters
        ----------
        wavelet:
            haar family: haar
            db family: db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12, db13, db14, db15, db16, db17, db18, db19, db20, db21, db22, db23, db24, db25, db26, db27, db28, db29, db30, db31, db32, db33, db34, db35, db36, db37, db38
            sym family: sym2, sym3, sym4, sym5, sym6, sym7, sym8, sym9, sym10, sym11, sym12, sym13, sym14, sym15, sym16, sym17, sym18, sym19, sym20
            coif family: coif1, coif2, coif3, coif4, coif5, coif6, coif7, coif8, coif9, coif10, coif11, coif12, coif13, coif14, coif15, coif16, coif17
            bior family: bior1.1, bior1.3, bior1.5, bior2.2, bior2.4, bior2.6, bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4, bior5.5, bior6.8
            rbio family: rbio1.1, rbio1.3, rbio1.5, rbio2.2, rbio2.4, rbio2.6, rbio2.8, rbio3.1, rbio3.3, rbio3.5, rbio3.7, rbio3.9, rbio4.4, rbio5.5, rbio6.8
            dmey family: dmey
            gaus family: gaus1, gaus2, gaus3, gaus4, gaus5, gaus6, gaus7, gaus8
            mexh family: mexh
            morl family: morl
            cgau family: cgau1, cgau2, cgau3, cgau4, cgau5, cgau6, cgau7, cgau8
            shan family: shan
            fbsp family: fbsp
            cmor family: cmor
        threshold:
            [0,1]
            higher removes more high signals
        level
        temporal_processing
        """

        super().__init__(temporal_processing)
        self.wavelet = wavelet
        self.threshold = threshold
        self.level = level

        try:
            import pywt
        except ImportError:
            raise ImportError("Please install pywt with `pip install pywt` to use Wavlet")

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        import pywt

        coeffs = pywt.wavedec(y, self.wavelet, level=self.level)

        # Apply soft thresholding to the detail coefficients
        threshold = np.sqrt(2 * np.log(len(y))) * self.threshold
        denoised_coeffs = [pywt.threshold(c, threshold, mode='soft') for c in coeffs]

        # Reconstruct the denoised signal
        denoised_signal = pywt.waverec(denoised_coeffs, wavelet=self.wavelet)

        return x, denoised_signal #[:-1]

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


def local_run():
    import matplotlib.pyplot as plt

    # Generate a noisy signal
    np.random.seed(42)
    time = np.linspace(0, 1, 500)
    signal = np.sin(2 * np.pi * 10 * time)  # Clean sine wave
    noise = np.random.normal(0, 0.5, time.shape)
    noisy_signal = signal + noise

    # Perform wavelet transform
    smoother = Wavelet(wavelet='sym4', level=4, threshold=0.2)
    time, denoised_signal = smoother.run(time, noisy_signal)

    # Plot the results
    plt.figure(figsize=(12, 6))
    plt.plot(time, noisy_signal, label='Noisy Signal', alpha=0.6)
    plt.plot(time, signal, label='Original Signal', linestyle='--')
    plt.plot(time, denoised_signal, label='Denoised Signal', linewidth=2)
    plt.legend()
    plt.title('Wavelet Denoising Example')
    plt.show()


if __name__ == "__main__":
    local_run()
