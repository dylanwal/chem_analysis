import numpy as np
import matplotlib.pyplot as plt
from skimage import data, restoration, util

import numpy as np
import matplotlib.pyplot as plt


from scipy.ndimage import minimum_filter, maximum_filter, grey_opening

from chem_analysis.processing.weigths.weights import DataWeight


def morphological_erosion(y: np.ndarray, window_size: int = 10):
    return minimum_filter(y, size=window_size, mode='nearest')


def morphological_dilation(y: np.ndarray, window_size: int = 10):
    return maximum_filter(y, size=window_size, mode='nearest')


def morphological_opening(y: np.ndarray, window_size: int = 10):
    return morphological_erosion(morphological_erosion(y), window_size)


def average_opening(y: np.ndarray, window_size: int = 10):
    return (morphological_dilation(grey_opening(y, window_size), window_size) +
            morphological_erosion(grey_opening(y, window_size), window_size)) / 2


# Example Data: Simulated Chromatogram with Baseline Drift
x = np.linspace(0, 100, 500)
signal = np.exp(-0.1 * (x - 50) ** 2) + np.exp(-2 * (x - 20) ** 2)+ 0.05 * x + 0.2 * np.sin(0.1 * x)  # Peak + Drift
noise = 0.02 * np.random.randn(len(x))  # Add noise
raw_signal = signal + noise

# Apply Adaptive Rolling Ball Algorithm with Correction
window = 45
# baseline = grey_opening(raw_signal, 45)
baseline = np.minimum(grey_opening(raw_signal, 45), average_opening(raw_signal, window))

# Plot Results
plt.figure(figsize=(10, 5))
plt.plot(x, raw_signal, label="Raw Signal", color='gray')
plt.plot(x, baseline, label=f"Estimated Baseline (r=)", linestyle="--", color="red")
plt.plot(x, raw_signal-baseline, label="Corrected Signal", color="blue")
plt.plot(x, average_opening(raw_signal, window), label="Corrected Signal", color="green")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Bidirectional Rolling Ball Baseline Correction")
plt.show()

# print(f"Optimized rolling ball radius: {final_radius}")
