import abc
import logging
from typing import Sequence

import numpy as np

import chem_analysis.utils.math as utils_math
from chem_analysis.library.chemical import Chemical
from chem_analysis.library.attribute import MassSpectrum
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.analysis.ms_analysis.ms_library import MSLibrary
import chem_analysis.utils.vector_similarity as vector_similarity

logger = logging.getLogger(__name__)


class MSScorer(abc.ABC):
    @abc.abstractmethod
    def __call__(self, mz: np.ndarray, library: np.ndarray, ms: np.ndarray) -> np.ndarray:
        ...


class ScorerMultiple(MSScorer):
    def __init__(self, scorers: list[MSScorer], mode: str = 'max'):
        self.scorers = scorers
        self.mode = mode

    def __call__(self, mz: np.ndarray, library: np.ndarray, ms: np.ndarray) -> np.ndarray:
        scores = []
        for scorer in self.scorers:
            scores.append(scorer(mz, library, ms))

        if self.mode == 'max':
            return np.max(scores, axis=0)
        if self.mode == 'min':
            return np.min(scores, axis=0)
        if self.mode == 'mean':
            return np.mean(scores, axis=0)
        if self.mode == 'median':
            return np.median(scores, axis=0)

        raise ValueError(f'Unknown mode {self.mode}')


class ScorerDot(MSScorer):
    def __init__(self, normalize: bool = True):
        """
        Larger Dot Product values indicate greater similarity.
        Dot Product is influenced by the length.
        If normalized; it is Cosine similarity.
        Cosine Similarity only considers the angle between vectors, regardless of their length/ magnitude.

            -1 = opposite directions
            0 = orthogonal directions
            1 = the identical
        """
        self.normalize = normalize

    def __call__(self, mz: np.ndarray, library: np.ndarray, ms: np.ndarray) -> np.ndarray:
        return vector_similarity.dot_cosine_similarity(library, ms, normalize=self.normalize)


class ScorerEarthMover(MSScorer):
    def __init__(self, normalize: bool = True):
        self.normalize = normalize

    def __call__(self, mz: np.ndarray, library: np.ndarray, ms: np.ndarray) -> np.ndarray:
        return vector_similarity.earth_movers_distance(mz, mz, library, ms, normalize=self.normalize)


class ScorerEuclidDistance(MSScorer):
    def __init__(self, normalize: bool = True):
        self.normalize = normalize

    def __call__(self, mz: np.ndarray, library: np.ndarray, ms: np.ndarray) -> np.ndarray:
        return vector_similarity.euclidean_distance(library, ms, self.normalize)


class MSFilter(abc.ABC):

    @abc.abstractmethod
    def __call__(self, chemicals: Sequence[Chemical], scores: np.ndarray) -> np.ndarray:
        ...


class FilterTopNMatches(MSFilter):
    def __init__(self, n: int):
        """ set all scores that are not in the top 'n' to zero. """
        self.n = n

    def __call__(self, chemicals: Sequence[Chemical], scores: np.ndarray) -> np.ndarray:
        sorted_score_index = np.argsort(-scores)
        scores[sorted_score_index[self.n:]] = 0
        return scores


class FilterMinScore(MSFilter):
    def __init__(self, score: int | float):
        """ set all scores less than provided score to zero. """
        self.score = score

    def __call__(self, chemicals: Sequence[Chemical], scores: np.ndarray) -> np.ndarray:
        mask = scores <= self.score
        scores[mask] = 0
        return scores


def search_by_ms(
        library: MSLibrary,
        ms: MSSignal,
        scorer_ms: MSScorer = ScorerDot,
        filter_ms: MSFilter | Sequence[MSFilter] = FilterMinScore(1),
) -> list:
    if not isinstance(filter_ms, list):
        filter_ms = [filter_ms]

    if len(library.ms_mass) == len(ms.x) and np.equal(ms.x, ms.y).all():
        ms_ = ms.y
    else:
        ms_ = utils_math.map_discrete_x_axis(library.ms_mass, ms.x, ms.y)

    scores = scorer_ms(library.ms_mass, library.ms_intensity, ms_)

    for filter_ in filter_ms:
        scores = filter_(library.chemicals, scores)

    labels = []
    sorted_score_index = np.argsort(-scores)
    for index in sorted_score_index:
        if scores[index] > 0:
            labels.append(library.chemicals[index])
    return labels


def search_by_ms_chemicals(
        chemicals: Sequence[Chemical],
        ms: MSSignal,
        scorer_ms: MSScorer = ScorerDot,
        filter_ms: MSFilter | Sequence[MSFilter] = FilterMinScore(1),
) -> list:
    if not isinstance(filter_ms, list):
        filter_ms = [filter_ms]

    scores = []
    for chemical in chemicals:
        for attr in chemical.attributes:
            if isinstance(attr, MassSpectrum):
                chem_ms = attr.value
                break
        else:
            logger.info(f'No mass spec found for {chemical}')
            scores.append(0)
            continue

        if len(chem_ms[:, 0]) == len(ms.x) and np.equal(ms.x, ms.y).all():
            mz, ms_1, ms_2 = ms.x, chem_ms[:, 1], ms.y
        else:
            mz, ms_1, ms_2 = map_discrete_x_axis(chem_ms[:, 0], chem_ms[:, 1], ms.x, ms.y)

        scores.append(scorer_ms(mz, ms_1, ms_2))

    if len(scores) == 0:
        return []

    scores = np.array(scores)

    for filter_ in filter_ms:
        scores = filter_(chemicals, scores)

    labels = []
    sorted_score_index = np.argsort(-scores)
    for index in sorted_score_index:
        if scores[index] > 0:
            labels.append(chemicals[index])
    return labels


def map_discrete_x_axis(
        x1: np.ndarray,
        y1: np.ndarray,
        x2: np.ndarray,
        y2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.arange(max(np.max(x1), np.max(x2))+1, dtype=int)
    y1_new = np.zeros_like(x)
    y2_new = np.zeros_like(x)

    y1_new[x1] = y1
    y2_new[x2] = y2

    return x, y1_new, y2_new
