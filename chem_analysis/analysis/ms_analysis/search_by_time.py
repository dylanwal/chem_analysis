
import numpy as np

from chem_analysis.analysis.peak import PeakContinuous
from chem_analysis.analysis.peak_result import ResultPeaks
from chem_analysis.analysis.ms_analysis.result_search import ResultCompoundSearch, PeakCompound
from chem_analysis.analysis.ms_analysis.ms_library import PickingLibrary


def get_top_n_matches(distance: np.ndarray, n: int = 1) -> np.ndarray:
    indices = np.argpartition(distance, n)[:n]
    return indices[np.argsort(-distance[indices])]


def search_by_retention_time(
        picking_library: PickingLibrary,
        peak_result: ResultPeaks,
        a_tolerance: int | float = 0.1,
        number_of_matches: int = 1,
) -> ResultCompoundSearch:
    peaks = []
    for peak in peak_result:
        peaks.append(search_by_retention_time_single(picking_library, peak, a_tolerance, number_of_matches))

    return ResultCompoundSearch(peaks)


def search_by_retention_time_single(
        picking_library: PickingLibrary,
        peak: PeakContinuous,
        a_tolerance: int | float = 0.1,
        number_of_matches: int = 1,
) -> PeakCompound:
    if not isinstance(peak, PeakContinuous):
        raise ValueError(f"Invalid 'Peak' type.\n peak type received: {type(peak)}")

    distance = np.abs(picking_library.retention_times - peak.max_x)
    index = get_top_n_matches(distance, number_of_matches)
    index = index[distance[index] < a_tolerance]

    if len(index) == 1:
        return PeakCompound(peak, picking_library.compounds[index[0]], distance[index[0]])
    elif len(index) == 0:
        return PeakCompound(peak, None, None)

    compounds = [picking_library.compounds[i] for i in index]
    distances = distance[index]
    return PeakCompound(peak, compounds, distances)


