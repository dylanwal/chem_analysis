import numpy as np

import chem_analysis as ca

rng = np.random.default_rng(seed=0)

x = np.linspace(0, 100, 500)
signal = np.exp(-0.1 * (x - 50) ** 2) + np.exp(-2 * (x - 20) ** 2) + 0.05 * x + 0.2 * np.sin(0.1 * x) + np.exp(-0.1 * (x - 80) ** 2)  # Peak + Drift
noise = 0.2 * rng.random(len(x))  # Add noise
y = signal + noise


def find_peaks_derivative(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    smoother = ca.p.smoothing.Gaussian(std=10)
    _, y_smooth = smoother.run(x, y)
    """Find peaks using the first and second derivatives."""
    dy = np.gradient(y_smooth, x)  # First derivative
    ddy = np.gradient(dy, x)  # Second derivative

    zero_crossings = np.where(np.diff(np.sign(dy)) < 0)[0]  # First derivative crosses zero downward
    peaks = [i for i in zero_crossings if ddy[i] < 0]

    return np.array(peaks, dtype=int)  # Peak x and y values


baseline = ca.p.baseline.MorphologicalAverage()
x, y = baseline.run(x, y)

x_peak, y_peak = find_peaks_derivative(x, y)

# Plot Results
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 5))
plt.plot(x, y, label="Raw Signal", color='gray')
plt.plot(x_peak, y_peak, 'o', label="new method", color="red")
# plt.plot(x, smoothed2, label="new method", linestyle="--", color="blue")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Bidirectional Rolling Ball Baseline Correction")
plt.show()