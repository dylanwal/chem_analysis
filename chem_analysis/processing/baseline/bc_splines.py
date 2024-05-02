import numpy as np
from scipy.interpolate import UnivariateSpline

from chem_analysis.processing.processing_method import Translation


class Spline(Translation):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 multiplier: float = 1,
                non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing)
        self.y_sub = y
        self.x_sub = x
        self.multiplier = multiplier

        # save results
        self.save_result = save_result
        self.baseline = None
        self.x = None
        self.data = None

    def run(self, x: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        baseline = self.get_baseline(x, data)
        data = data - baseline

        if self.save_result:
            self.baseline = baseline
            self.x = x
            self.data = data

        return x, data

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # if len(self.y_sub) == len(y):
            # return self.multiplier * self.y_sub

        # spline = UnivariateSpline(x, y, s=0)

        raise NotImplementedError()  # TODO: x-interpolation

    def get_baseline2D(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)

        for i in range(z.shape[0]):
            baseline[i, :] = self.get_baseline(x, z[i, :])

        return baseline
