from typing import Iterable

import numpy as np

from chem_analysis.analysis.peak import PeakContinuous
from chem_analysis.analysis.peak_picking.picking_peak import PeakForPicking
from chem_analysis.analysis.mass_spec_analysis.criteria import Criteria


class PickingLibrary:
    def __init__(self,
                 peaks: Iterable[PeakForPicking],
                 criteria: Criteria | Iterable[Criteria] = None,
                 all_criteria: bool = False
                 ):
        self.peaks = peaks
        self.criteria = criteria
        self.all_criteria = all_criteria

    def __iter__(self):
        return iter(self.peaks)

    def search_lib(self, signal_peak: PeakContinuous) -> PeakForPicking | None:
        for peak in self.peaks:
            possible_matches = []
            if self.within_tolerance(peak, signal_peak):
                possible_matches.append(peak)
            if possible_matches:
                diff = [abs(signal_peak.max_loc - possible_peak.pos_x) for possible_peak in possible_matches]
                return possible_matches[np.argmin(diff)]

        return None

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
