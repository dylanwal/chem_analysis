import numpy as np

from chem_analysis.processing.processing_method import Resampling


class ResampleDtypeNormalize(Resampling):

    def __init__(self,
                 dtype: np.dtype = None,
                 ):
        """

        Parameters
        ----------
        dtype:

        """
        super().__init__()
        self.dtype = dtype

    @property
    def dtype_max(self):
        return np.iinfo(self.dtype).max

    @property
    def dtype_min(self):
        return np.iinfo(self.dtype).min

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x, (y / np.max(y) * self.dtype_max).astype(self.dtype)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return x, y, (z / np.max(z) * self.dtype_max).astype(self.dtype)

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return x, y, z, (w / np.max(w) * self.dtype_max).astype(self.dtype)
