

class LocalMax(DiscoveryMethods):
    def __init__(self):
        pass

    def run(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        from scipy.signal._peak_finding_utils import _local_maxima_1d
        index, _, _ = _local_maxima_1d(y)
        return index


class Derivative(DiscoveryMethods):
    def __init__(self):
        pass

    def run(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Find peaks using the first and second derivatives."""
        dy = np.gradient(y, x)  # First derivative
        ddy = np.gradient(dy, x)  # Second derivative

        zero_crossings = np.where(np.diff(np.sign(dy)) < 0)[0]  # First derivative crosses zero downward
        peaks = [i for i in zero_crossings if ddy[i] < 0]

        return np.array(peaks, dtype=int)  # Peak x and y values