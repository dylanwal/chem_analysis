from typing import Iterable

from chem_analysis.analysis.peak_picking.library_search.criteria import Criteria, CriteriaAbsoluteRangeAnd


class PeakForPicking:
    DEFAULT_CRITERIA = CriteriaAbsoluteRangeAnd(atol_x=(0.01, 0.01))

    def __init__(self,
                 pos_x: int | float,
                 pos_y: int | float = None,
                 criteria: Criteria | Iterable[Criteria] = None,
                 all_criteria: bool = False
                 ):
        """

        Parameters
        ----------
        pos_x
        pos_y
        criteria
        all_criteria:
            True: all must criteria must be true to include
            False: one criteria producing True is sufficient
        """
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.criteria = criteria
        self.all_criteria = all_criteria

    def within_tolerance(self, pos_x: int | float, pos_y: int | float = None) -> bool:
        if self.criteria is None:
            criteria = self.DEFAULT_CRITERIA
        else:
            criteria = self.criteria

        if not isinstance(criteria, Iterable):
            criteria = [criteria]

        for criteria in criteria:
            result = criteria.evaluate(self, pos_x, pos_y)
            if self.all_criteria and result:
                continue
            if result:
                return True
        return False


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

    def __contains__(self, item):
        for peak in self.peaks:
            if self.within_tolerance(peak, item):
                return True
        return False

    def within_tolerance(self, lib_peak: PeakForPicking, signal_peak: Peak) -> bool:
        if self.criteria is None:
            return True

        criteria = self.criteria
        if not isinstance(criteria, Iterable):
            criteria = [criteria]

        for criteria in criteria:
            result = criteria.evaluate(lib_peak, signal_peak)
            if self.all_criteria and result:
                continue
            if result:
                return True

        peak.within_tolerance(pos_x, pos_y)

        return False
