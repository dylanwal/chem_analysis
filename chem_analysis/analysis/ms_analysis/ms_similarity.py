"""
Relevant paper: DOI: 10.1016/1044-0305(94)87009-8
-> dot prodcut outperforms Euclidean distance and probability-matching


"""

import numpy as np

from chem_analysis.mass_spec.ms_signal import MSSignal
import chem_analysis.utils.math as utils_math


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def normalize_matrix(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.where(norms != 0, norms, 1)


def dot_product_signals(signal1: MSSignal, signal2: MSSignal) -> float:
    """

    Parameters
    ----------
    signal1
    signal2

    Returns
    -------
    values: [0, 1]
        0 = not similar
        1 = the identical (note identical vectors may not return 1 due to numerical roundoff errors).
    """
    if np.all(np.isclose(signal1.x, signal2.x, rtol=0.01)):
        y_1, y_2 = signal1.y, signal2.y
    else:
        # unify mz
        max_mz = int(max(signal1.mz.max(), signal2.mz.max()))
        x = np.arange(max_mz, dtype=utils_math.min_int_dtype(max_mz, 0))
        y_1 = utils_math.map_discrete_x_axis(x, np.round(signal1.mz), signal1.y)
        y_2 = utils_math.map_discrete_x_axis(x, np.round(signal2.mz), signal2.y)

    y_1 = normalize_vector(y_1)
    y_2 = normalize_vector(y_2)
    return np.dot(y_1, y_2)


def dot_product_matrix(signals: np.ndarray, target_ms: np.ndarray) -> np.ndarray:
    """

    Parameters
    ----------
    signals: [x, y]
        2d matrix of signals that you want compare
    target_ms: [x]

    Returns
    -------
    [y]
    """
    return np.dot(normalize_matrix(signals), normalize_vector(target_ms))
