import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


def find_derivative_bounds(x, y, cutoff: float = 0.1, filter =None):
    # Find the peak
    peaks, _ = find_peaks(y)

    peak_index = peaks[np.argmax(y[peaks])]

    # Calculate the first derivative
    dy = np.gradient(y, x)

    # Find left bound
    max_ = 0
    left_index =None
    for i in range(peak_index):
        dy_ = dy[peak_index-i]
        if dy_ > max_:
            max_ = dy_
        else:
            if max_*cutoff > dy_:
                left_index = peak_index-i
                break

    if left_index is None:
        left_bound = x[0]
    else:
        left_bound = x[left_index]

    # Find right bound
    min_ = 0
    right_index = None
    for i in range(len(dy) - peak_index):
        dy_ = dy[peak_index+i+1]
        if dy_ < min_:
            min_ = dy_
        else:
            if min_*cutoff < dy_:
                right_index = peak_index+i
                break

    if right_index is None:
        right_bound = x[-1]
    else:
        right_bound = x[right_index]

    return left_bound, right_bound


def run_local():
    # Example data
    x = np.linspace(0, 10, 1000)
    y = np.exp(-(x - 5)**2) + 0.1 * np.random.random(x.size)

    # Find bounds using the derivative method
    left_bound, right_bound = find_derivative_bounds(x, y)
    print(f"Bounds of integration (derivative method): [{left_bound}, {right_bound}]")
    # Plot the original data and the fitted Gaussian
    plt.plot(x, y, label='Original Data')
    plt.axvline(x=left_bound, color='r', linestyle=':', label='Left Bound')
    plt.axvline(x=right_bound, color='g', linestyle=':', label='Right Bound')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run_local()
