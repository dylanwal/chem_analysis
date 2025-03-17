from typing import Protocol

import numpy as np
from scipy.stats import entropy, wasserstein_distance


class SimilarityFunction(Protocol):
    def __call__(self, vector1: np.ndarray, vector2: np.ndarray, *args, **kwargs) -> float | np.ndarray:
        ...


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector, axis=-1, keepdims=True)  # Keep dimensions consistent
    return vector / norm  # Safe broadcasting


def dot_cosine_similarity(vec1: np.ndarray, vec2: np.ndarray, normalize: bool = True) -> float | np.ndarray:
    """
    Larger Dot Product values indicate greater similarity.
        Dot Product is influenced by the length.
    If normalized; it is Cosine similarity.
        Cosine Similarity only considers the angle between vectors, regardless of their length/ magnitude.

    -1 = opposite directions
    0 = orthogonal directions
    1 = the identical

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray
    normalize: bool
        make vectors unit vectors
        output limited to [-1, 1]

    Returns
    -------
    value: float
        [-inf, inf] or [-1, 1] when normalize=True

    """
    if normalize:
        vec1 = normalize_vector(vec1)
        vec2 = normalize_vector(vec2)
    return np.dot(vec1, vec2)


def earth_movers_distance(mz1: np.ndarray, mz2: np.ndarray, vec1: np.ndarray, vec2: np.ndarray, normalize: bool = True) -> float |np.ndarray:
    """
    Compute the Earth Mover's Distance (EMD) between two mass spectra.

    Parameters:
    - mz1, intensity1: Arrays representing the m/z values and intensities of spectrum 1.
    - mz2, intensity2: Arrays representing the m/z values and intensities of spectrum 2.

    Returns:
    - EMD distance between the two spectra.
    """
    if normalize:
        vec1 = normalize_vector(vec1)
        vec2 = normalize_vector(vec2)
    return np.array([wasserstein_distance(mz1, mz2, u_weights=row, v_weights=vec2) for row in vec1])


def euclidean_distance(vec1: np.ndarray, vec2: np.ndarray, normalize: bool = True) -> float | np.ndarray:
    """
    Measures the straight-line (L2) distance between two vectors.

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray
    normalize: bool

    Returns
    -------
    value: float

    """
    if normalize:
        vec1 = normalize_vector(vec1)
        vec2 = normalize_vector(vec2)
    return np.sqrt(np.sum((vec1 - vec2) ** 2, axis=1))


def manhattan_distance(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    L1 distance or "taxicab" distance, it sums the absolute differences between corresponding elements

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray

    Returns
    -------
    value: float

    """
    return np.sum(np.abs(vec1 - vec2))


def jaccard_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Measures the size of the intersection divided by the size of the union of the two sets, suitable for binary vectors or sets.
    Useful in cases where only the presence/absence of elements matters (e.g., document similarity).

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray

    Returns
    -------
    value: float

    """
    intersection = np.minimum(vec1, vec2).sum()
    union = np.maximum(vec1, vec2).sum()
    return intersection / union if union != 0 else 0


def pearson_correlation(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Measures linear correlation between two vectors, scaled between -1 and 1

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray

    Returns
    -------
    value: float

    """
    return np.corrcoef(vec1, vec2)[0, 1]


def hamming_distance(vec1: np.ndarray, vec2: np.ndarray, binary: bool = True) -> float:
    """
    Counts the number of positions at which the corresponding elements differ, mainly for binary or categorical vectors.
    Useful in genetic sequences, error detection/correction, and binary data.

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray
    binary: bool
        turns values into 0 if zero or 1 is non-zero before computing hamming distance

    Returns
    -------
    value: float
    """
    if binary:
        vec1 = np.abs(vec1) > 0
        vec2 = np.abs(vec2) > 0
        return np.sum(np.logical_and(vec1, vec2))

    return np.sum(vec1 != vec2) / len(vec1)


def mahalanobis_distance(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Measures distance between two points while considering correlations within the data.
    Useful for identifying multivariate outliers and data with correlated features.

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray

    Returns
    -------
    value: float
    """
    cov_matrix = np.cov(np.vstack([vec1, vec2]).T)
    diff = vec1 - vec2
    return np.sqrt(np.dot(np.dot(diff.T, np.linalg.inv(cov_matrix)), diff))


def kl_divergence(vec1: np.ndarray, vec2: np.ndarray, offset: float = 1) -> float:
    """
    Measures how one probability distribution diverges from a second, expected distribution.
    Often used in probability and information theory for distributions.

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray
    offset: float
        used when zero is present in denominator
        scalar for min_value

    Returns
    -------
    value: float
    """
    vec1 = np.asarray(vec1, dtype=np.float64)
    if np.count_nonzero(vec2) > 0:
        # avoid dividing by zero
        min_value = np.min(vec2[np.abs(vec2) > 0])
        vec2 = np.asarray(vec2, dtype=np.float64) + min_value * offset
    else:
        vec2 = np.asarray(vec2, dtype=np.float64)
    vec1 /= np.sum(vec1)
    vec2 /= np.sum(vec2)
    return entropy(vec1, vec2)


def dice_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Dice similarity (or Dice coefficient) is a metric similar to Jaccard similarity but is often more sensitive in
    applications like comparing binary or sparse vectors.
    It's especially useful in comparing sets or binary vectors where elements are either "present" (1) or "absent" (0).

    Parameters
    ----------
    vec1: np.ndarray
    vec2: np.ndarray

    Returns
    -------
    value: float
    """
    intersection = np.sum(vec1 * vec2)
    return (2 * intersection) / (np.sum(vec1) + np.sum(vec2)) if (np.sum(vec1) + np.sum(vec2)) != 0 else 0
