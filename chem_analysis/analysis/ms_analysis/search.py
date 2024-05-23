from typing import Iterable

from chem_analysis.analysis.picking_result import ResultPeaks
from chem_analysis.analysis.ms_analysis.picking_library import PickingLibrary
from chem_analysis.analysis.ms_analysis.search_result import CompoundSearchResult
from chem_analysis.analysis.ms_analysis.criteria import Criteria, CriteriaAbsoluteRangeAnd


DEFAULT_CRITERIA = CriteriaAbsoluteRangeAnd(atol_x=(0.1, 0.1))


def search_by_retention_time(
        picking_library: PickingLibrary,
        peak_result: ResultPeaks,
        criteria: Criteria | Iterable[Criteria] = None,
        number_of_matches: int = 1,
) -> CompoundSearchResult:
    for peak in peak_result:
        possible_matches = []
        if self.within_tolerance(peak, signal_peak):
            possible_matches.append(peak)
        if possible_matches:
            diff = [abs(signal_peak.max_loc - possible_peak.pos_x) for possible_peak in possible_matches]
            return possible_matches[np.argmin(diff)]

    return CompoundSearchResult()


def search_by_retention_time_peak(picking_library: PickingLibrary, criteria: Criteria, peak) -> CompoundSearchResult:
    picking_library.retention_times


def within_tolerance(self, lib_peak: PeakForPicking, signal_peak: PeakContinuous) -> bool:
    if self.criteria is not None:
        criteria = self.criteria
        if not isinstance(criteria, Iterable):
            criteria = [criteria]

        for criteria_ in criteria:
            result = criteria_.evaluate(lib_peak, signal_peak)
            if self.all_criteria and result:
                continue
            if result:
                return True

    return lib_peak.within_tolerance(signal_peak)

# def get_n_nearest_compounds(self, retention_time: float, n: int = 1) -> list[Compound]:
#     distance = []
#     for i, peak in enumerate(self.peaks):
#         distance[i] = abs(peak.retention_time - retention_time)
#
#     sort_index = np.argsort(distance)
#     compounds = []
#     for i in range(n):
#         compounds.append(self.compounds[sort_index[i]])
#     return compounds