"""
Relevant paper: DOI: 10.1016/1044-0305(94)87009-8
-> dot product outperforms Euclidean distance and probability-matching

"""

from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.mass_spec.ms_utils import unify_ms_signals
import chem_analysis.utils.vector_similarity as vector_similarity


def ms_similarity(
        signal1: MSSignal,
        signal2: MSSignal,
        method: vector_similarity.SimilarityFunction = vector_similarity.dot_cosine_similarity
) -> float:
    """

    Parameters
    ----------
    signal1
    signal2
    method: vector_similarity.SimilarityFunction
        method used to compare two MS signals

    Returns
    -------
    values: float
    """
    ms1, ms2 = unify_ms_signals(signal1, signal2)
    return method(ms1.y, ms2.y)

