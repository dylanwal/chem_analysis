from typing import Protocol

import numpy as np
from scipy.stats import entropy, wasserstein_distance


class SimilarityFunction(Protocol):
    def __call__(self, vector1: np.ndarray, vector2: np.ndarray, *args, **kwargs) -> float | np.ndarray:
        ...


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector, axis=-1, keepdims=True)  # Keep dimensions consistent
    return vector / norm  # Safe broadcasting


# def dot_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float | np.ndarray:
#     """
#     Larger Dot Product values indicate greater similarity.
#         Dot Product is influenced by the length.
#
#     Parameters
#     ----------
#     vec1: np.ndarray [n, m] or [n]
#         This should be the library of values you want to compare.
#     vec2: np.ndarray [n]
#         The vector you want to compare.
#
#     Returns
#     -------
#     distance:
#         [-inf, inf]
#         neg = opposite directions
#         0 = orthogonal directions
#         pos = the identical
#
#     """
#     return np.dot(vec1, vec2) # not a distance


def cosine_distance(vec1: np.ndarray, vec2: np.ndarray) -> float | np.ndarray:
    """
    Cosine distance
        Only considers the angle between vectors, regardless of their length/ magnitude.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance:
        [0, 2] when normalize=True
        2 = opposite directions
        1 = orthogonal directions
        0 = the identical
    """
    vec1 = normalize_vector(vec1)
    vec2 = normalize_vector(vec2)
    return 1 - np.dot(vec1, vec2)


def earth_movers_distance(
        vec1: np.ndarray,
        vec2: np.ndarray,
        mz1: np.ndarray = None,
        mz2: np.ndarray = None
) -> float | np.ndarray:
    """
    Compute the Earth Mover's Distance (EMD) between two mass spectra.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.
    mz1: np.ndarray [n]
    mz2: np.ndarray [n]

    Returns
    -------
    distance:
        small values means vectors are similar

    """
    if len(vec1.shape) == 1:
        if mz1 is None or mz2 is None:
            mz1 = mz2 = np.arange(vec1.shape[0])
        return wasserstein_distance(mz1, mz2, u_weights=vec1, v_weights=vec2)

    if mz1 is None or mz2 is None:
        mz1 = mz2 = np.arange(vec1.shape[1])
    if len(vec1.shape) == 1:
        return wasserstein_distance(mz1, mz2, u_weights=vec1, v_weights=vec2)
    return np.array([wasserstein_distance(mz1, mz2, u_weights=row, v_weights=vec2) for row in vec1])


def euclidean_distance(vec1: np.ndarray, vec2: np.ndarray, normalize: bool = True) -> float | np.ndarray:
    """
    Measures the straight-line (L2) distance between two vectors.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.
    normalize: bool

    Returns
    -------
    distance:
        small values means vectors are similar

    """
    if normalize:
        vec1 = normalize_vector(vec1)
        vec2 = normalize_vector(vec2)
    if len(vec1.shape) == 1:
        return np.sqrt(np.sum((vec1 - vec2) ** 2))
    return np.sqrt(np.sum((vec1 - vec2) ** 2, axis=1))


def manhattan_distance(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    L1 distance or "taxicab" distance, it sums the absolute differences between corresponding elements

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance:
        small values means vectors are similar

    """
    if len(vec1.shape) == 1:
        return np.sum(np.abs(vec1 - vec2))
    return np.sum(np.abs(vec1 - vec2), axis=1)


def jaccard_similarity(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray | float:
    """
    Measures the size of the intersection divided by the size of the union of the two sets, suitable for binary vectors or sets.
    Useful in cases where only the presence/absence of elements matters (e.g., document similarity).

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance:
        0: perfect match
        1: no match

    """
    vec1 = normalize_vector(vec1)
    vec2 = normalize_vector(vec2)

    if len(vec1.shape) == 1:
        intersection = np.sum(np.minimum(vec1, vec2))
        union = np.sum(np.maximum(vec1, vec2))
        return 1 - intersection / union if union != 0 else 1

    intersection = np.sum(np.minimum(vec1, vec2), axis=1)
    union = np.sum(np.maximum(vec1, vec2), axis=1)
    nonzero_index = np.nonzero(union)[0]
    if len(nonzero_index) == len(intersection):
        return 1 - intersection / union

    out = np.zeros_like(intersection)
    out[nonzero_index] = intersection[nonzero_index] / union[nonzero_index]
    return 1 - out


def pearson_correlation(vec1: np.ndarray, vec2: np.ndarray) -> float | np.ndarray:
    """
    Measures linear correlation between two vectors, scaled between -1 and 1

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance:

    """
    if len(vec1.shape) == 1:
        return 1 - np.corrcoef(vec1, vec2)[0, 1]
    return 1 - np.corrcoef(vec1, vec2)[-1, :-1]


def hamming_distance(vec1: np.ndarray, vec2: np.ndarray, binary: bool = True) -> float:
    """
    Counts the number of positions at which the corresponding elements differ, mainly for binary or categorical vectors.
    Useful in genetic sequences, error detection/correction, and binary data.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.
    binary: bool
        turns values into 0 if zero or 1 is non-zero before computing hamming distance

    Returns
    -------
    distance:
    """
    vec1 = np.abs(vec1) > 0
    vec2 = np.abs(vec2) > 0

    if len(vec1.shape) == 1:
        if binary:
            return np.sum(np.logical_and(vec1, vec2))

        return np.sum(vec1 != vec2) / len(vec1)

    if binary:
        return np.sum(np.logical_and(vec1, vec2), axis=1)

    return np.sum(vec1 != vec2, axis=1) / len(vec1)


def mahalanobis_distance(vec1: np.ndarray, vec2: np.ndarray) -> float | np.ndarray:
    """
    Measures distance between two points while considering correlations within the data.
    Useful for identifying multivariate outliers and data with correlated features.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance
    """
    vec1 = np.atleast_2d(vec1)
    vec2 = np.asarray(vec2).reshape(1, -1)

    if vec1.shape[1] != vec2.shape[1]:
        raise ValueError("vec1 and vec2 must have the same number of features")

    # Compute covariance matrix from vec1 only
    cov = np.cov(vec1, rowvar=False)
    cov_inv = np.linalg.pinv(cov)

    diffs = vec1 - vec2
    dists = np.sqrt(np.einsum("ij,jk,ik->i", diffs, cov_inv, diffs))
    return dists[0] if dists.shape[0] == 1 else dists

    # import scipy.distance as distance
    # out = np.empty(vec1.shape[0])
    # cov_ = np.cov(vec1, rowvar=False)
    # for i in range(vec1.shape[0]):
    #     out[i] = distance.mahalanobis(vec1[i], vec2.reshape(vec2.size), cov_)
    #
    # return out


def kl_divergence(vec1: np.ndarray, vec2: np.ndarray, offset: float = 1) -> float | np.ndarray:
    """
    Kullback-Leibler divergence
    Measures how one probability distribution diverges from a second, expected distribution.
    Often used in probability and information theory for distributions.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance
    """
    if np.count_nonzero(vec2) > 0:
        # avoid dividing by zero
        min_value = np.min(vec2[np.abs(vec2) > 0])
        vec2 = np.asarray(vec2, dtype=np.float64) + min_value * offset
    else:
        vec2 = np.asarray(vec2, dtype=np.float64)
    vec2 = vec2/np.sum(vec2)

    if len(vec1.shape) == 1:
        vec1 = vec1/np.sum(vec1)
        return entropy(vec1, vec2)

    vec1 = vec1/np.sum(vec1, axis=1, keepdims=True)
    return np.array(tuple(entropy(vec, vec2) for vec in vec1))


def dice_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float | np.ndarray:
    """
    Dice similarity (or Dice coefficient) is a metric similar to Jaccard similarity but is often more sensitive in
    applications like comparing binary or sparse vectors.
    It's especially useful in comparing sets or binary vectors where elements are either "present" (1) or "absent" (0).

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance
    """
    vec1 = normalize_vector(vec1)
    vec2 = normalize_vector(vec2)

    if len(vec1.shape) == 1:
        intersection = np.sum(vec1 * vec2)
        return 1 - (2 * intersection) / (np.sum(vec1) + np.sum(vec2)) if (np.sum(vec1) + np.sum(vec2)) != 0 else 0

    intersection = np.sum(vec1 * vec2, axis=1)
    nonzero_index = np.nonzero(np.sum(vec1, axis=1) + np.sum(vec2))[0]
    if len(nonzero_index) == len(intersection):
        return 1 - (2 * intersection) / (np.sum(vec1, axis=1) + np.sum(vec2))

    out = np.zeros_like(intersection)
    out[nonzero_index] = (2 * intersection[nonzero_index]) / (np.sum(vec1, axis=1)[nonzero_index] + np.sum(vec2))
    return 1 - out


def weighted_recall_score(vec1: np.ndarray, vec2: np.ndarray, beta: float = 0.1) -> np.ndarray | float:
    """
    Computes a weighted recall-based score.
    This function measures how much of the expected signal is present while applying a soft penalty for extra elements.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.
    beta:
        Penalty factor for extra elements.

    Returns
    -------
    distance
    """
    def score(single_vec1, single_vec2):
        expected = set(single_vec1)
        observed = set(single_vec2)

        if not expected:
            return 1.0 if not observed else 0.0

        true_positive = len(expected & observed)
        false_positive = len(observed - expected)

        recall = true_positive / len(expected)
        penalty = beta * (false_positive / len(observed))

        return recall - np.minimum(penalty, recall)

    if vec1.ndim == 1:
        return 1 - score(vec1, vec2)

    return 1 - np.array([score(row, vec2) for row in vec1])


def dot_product_soft_norm(vec1: np.ndarray, vec2: np.ndarray, gamma: float = 0.5) -> float | np.ndarray:
    """
    Computes a dot product similarity with soft normalization.
    This function evaluates similarity while softly penalizing excess signal.

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.
    gamma:
        Penalty factor for extra signal

    Returns
    -------
    distance
    """
    vec1 = normalize_vector(vec1)
    vec2 = normalize_vector(vec2)

    vec1 = np.atleast_2d(vec1)
    vec2 = np.asarray(vec2).reshape(1, -1)

    if vec1.shape[1] != vec2.shape[1]:
        raise ValueError("vec1 and vec2 must have the same number of features")

    dot_products = np.sum(vec1 * vec2, axis=1)
    norm_expected = np.linalg.norm(vec1, axis=1)
    norm_observed = np.linalg.norm(vec2)

    with np.errstate(divide='ignore', invalid='ignore'):
        alignment_score = np.divide(dot_products, norm_expected, out=np.zeros_like(dot_products, dtype=float), where=norm_expected != 0)

        diff = vec1 - vec2
        penalty = np.linalg.norm(diff, axis=1)
        penalty = gamma * np.divide(penalty, norm_observed if norm_observed else 1)

    score = alignment_score - penalty
    return 1 - score[0] if score.shape[0] == 1 else 1 - score


def dtw_distance(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """
    Compute the DTW distance between one or more sequences (vec1) and a reference sequence (vec2).

    Parameters
    ----------
    vec1: np.ndarray [n, m] or [n]
        This should be the library of values you want to compare.
    vec2: np.ndarray [n]
        The vector you want to compare.

    Returns
    -------
    distance
    """
    vec1 = normalize_vector(vec1)
    vec2 = normalize_vector(vec2)

    def _dtw_1d(a: np.ndarray, b: np.ndarray) -> float:
        n, m = len(a), len(b)
        dtw_matrix = np.full((n + 1, m + 1), np.inf)
        dtw_matrix[0, 0] = 0.0

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                cost = abs(a[i - 1] - b[j - 1])
                dtw_matrix[i, j] = cost + min(
                    dtw_matrix[i - 1, j],    # insertion
                    dtw_matrix[i, j - 1],    # deletion
                    dtw_matrix[i - 1, j - 1] # match
                )

        return dtw_matrix[n, m]

    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2).flatten()

    if vec1.ndim == 1:
        return _dtw_1d(vec1, vec2)
    elif vec1.ndim == 2:
        return np.array([_dtw_1d(row, vec2) for row in vec1])
    else:
        raise ValueError("vec1 must be a 1D or 2D numpy array")
