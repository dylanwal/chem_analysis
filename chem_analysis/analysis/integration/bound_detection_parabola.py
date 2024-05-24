
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


def find_parabola_bounds(x, y, window_size=5):
    # Find the peak
    peaks, _ = find_peaks(y)
    if not peaks:
        return None, None

    peak_index = peaks[np.argmax(y[peaks])]

    # Initialize bounds
    left_bound = None
    right_bound = None

    # Check left side
    for i in range(peak_index, window_size - 1, -1):
        window_x = x[i - window_size + 1: i + 1]
        window_y = y[i - window_size + 1: i + 1]
        coeffs = np.polyfit(window_x, window_y, 2)
        second_derivative = 2 * coeffs[0]  # second derivative of ax^2 + bx + c is 2a
        if second_derivative > 0:
            left_bound = x[i]
            break

    # Check right side
    for i in range(peak_index, len(x) - window_size):
        window_x = x[i: i + window_size]
        window_y = y[i: i + window_size]
        coeffs = np.polyfit(window_x, window_y, 2)
        second_derivative = 2 * coeffs[0]
        if second_derivative > 0:
            right_bound = x[i]
            break

    return left_bound, right_bound

# Example data
x = np.linspace(0, 10, 1000)
y = np.exp(-(x - 5)**2)

# Find bounds using the parabola fitting method
left_bound, right_bound = find_parabola_bounds(x, y, window_size=30)
print(f"Bounds of integration (parabola fitting method): [{left_bound}, {right_bound}]")

# Plot the original data and the bounds
plt.plot(x, y, label='Original Data')
plt.axvline(x=left_bound, color='r', linestyle=':', label='Left Bound')
plt.axvline(x=right_bound, color='g', linestyle=':', label='Right Bound')
plt.legend()
plt.show()
