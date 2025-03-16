import abc

import numpy as np

import chem_analysis.utils.math as utils_math
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.analysis.ms_analysis.ms_library import MSLibrary
import chem_analysis.utils.vector_similarity as vector_similarity


class MSFilter(abc.ABC):
    def __init__(self):
        ...

    @abc.abstractmethod
    def run(self, lib: MSLibrary, scores: np.ndarray) -> np.ndarray:
        ...


class TopNMatches:
    def __init__(self, n: int):
        self.n = n

    def run(self, lib: MSLibrary, scores: np.ndarray) -> np.ndarray:
        indices = np.argpartition(scores, -self.n)[-self.n:]
        return indices[np.argsort(-scores[indices])]


def search_by_ms(
        library: MSLibrary,
        ms: MSSignal,
        scorer_ms: vector_similarity.SimilarityFunction = vector_similarity.dot_cosine_similarity,
        filter_ms: MSFilter = TopNMatches(n=1),
) -> list | list[list]:
    if len(library.ms_mass) == len(ms.x) and np.equal(ms.x, ms.y).all():
        ms_ = ms.y
    else:
        ms_ = utils_math.map_discrete_x_axis(library.ms_mass, ms.x, ms.y)

    scores = scorer_ms(library.ms_intensity, ms_)
    index = filter_ms.run(library, scores)
    labels = []
    for i in index:
        if index == -1:
            labels.append(None)
        else:
            if isinstance(i, list):
                labels.append(library.chemicals[ii] for ii in i)
            else:
                labels.append(library.chemicals[i])

    return labels
