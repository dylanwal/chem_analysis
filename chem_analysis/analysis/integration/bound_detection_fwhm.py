
def find_fwhm_bounds(x, y):
    # Find the peak
    peaks, _ = find_peaks(y)
    if not peaks:
        return None, None

    peak_index = peaks[np.argmax(y[peaks])]
    y_peak = y[peak_index]
    y_half = y_peak / 2

    # Find left bound
    left_index = np.where(y[:peak_index] <= y_half)[0]
    if left_index.size > 0:
        left_bound = x[left_index[-1]]
    else:
        left_bound = x[0]

    # Find right bound
    right_index = np.where(y[peak_index:] <= y_half)[0]
    if right_index.size > 0:
        right_bound = x[peak_index + right_index[0]]
    else:
        right_bound = x[-1]

    return left_bound, right_bound

# Example data
x = np.linspace(0, 10, 1000)
y = np.exp(-(x - 5)**2)

# Find bounds
left_bound, right_bound = find_fwhm_bounds(x, y)
print(f"Bounds of integration: [{left_bound}, {right_bound}]")
