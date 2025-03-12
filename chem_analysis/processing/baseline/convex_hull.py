from typing import Iterable

import numpy as np
import scipy.spatial

from chem_analysis.processing.processing_method import Baseline


def rubber_band(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    xy = np.column_stack((x, y))
    v = scipy.spatial.ConvexHull(xy).vertices
    v = np.roll(v, -v.argmin())
    v = v[:(v.argmax()+1)]

    return np.interp(x, x[v], y[v])


class ConvexHull(Baseline):
    def __init__(self,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return rubber_band(x, y)


def marching_lines(x: np.ndarray, y: np.ndarray, len: float = 1, max_angle: float = 0.2) -> np.ndarray:
    """

    Parameters
    ----------
    x
    y
    len
    max_angle

    Returns
    -------

    """
    #start on one end; start with line pointed down and rotate up till the line hits a point. then use the
    # end of the line as the  next point and repeat.
    # TODO: