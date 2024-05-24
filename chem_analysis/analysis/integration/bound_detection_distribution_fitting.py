

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define a Gaussian function
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def find_statistical_bounds(x, y):
    # Guess initial parameters for Gaussian fitting
    mean_guess = x[np.argmax(y)]
    amplitude_guess = np.max(y)
    stddev_guess = np.std(x)

    # Fit the Gaussian model
    popt, _ = curve_fit(gaussian, x, y, p0=[amplitude_guess, mean_guess, stddev_guess])

    amplitude, mean, stddev = popt

    # Calculate bounds as mean ± 3*stddev
    left_bound = mean - 3 * stddev
    right_bound = mean + 3 * stddev

    return left_bound, right_bound, amplitude, mean, stddev

# Example data
x = np.linspace(0, 10, 1000)
y = np.exp(-(x - 5)**2)

# Find bounds using the statistical method
left_bound, right_bound, amplitude, mean, stddev = find_statistical_bounds(x, y)
print(f"Bounds of integration (statistical method): [{left_bound}, {right_bound}]")

# Plot the original data and the fitted Gaussian
plt.plot(x, y, label='Original Data')
plt.plot(x, gaussian(x, amplitude, mean, stddev), label='Fitted Gaussian', linestyle='--')
plt.axvline(x=left_bound, color='r', linestyle=':', label='Left Bound')
plt.axvline(x=right_bound, color='g', linestyle=':', label='Right Bound')
plt.legend()
plt.show()
