import copy
from typing import Iterable
from collections import OrderedDict

from chem_analysis.analysis.peak import PeakContinuous, PeakBounded, PeakParent
from chem_analysis.analysis.peak_picking.library_search.criteria import Criteria, CriteriaAbsoluteRangeAnd


class PeakForPicking:
    DEFAULT_CRITERIA = CriteriaAbsoluteRangeAnd(atol_x=(0.1, 0.1))

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

    def within_tolerance(self, peak: PeakContinuous) -> bool:
        if self.criteria is None:
            criteria = self.DEFAULT_CRITERIA
        else:
            criteria = self.criteria

        if not isinstance(criteria, Iterable):
            criteria = [criteria]

        for criteria in criteria:
            result = criteria.evaluate(self, peak)
            if self.all_criteria and result:
                continue
            if result:
                return True
        return False

    def create_peak(self, peak: PeakContinuous):
        return copy.copy(peak)


from chem_analysis.gc_lc.gc_library import Compound


class PeakForPickingCompound(PeakForPicking):
    def __init__(self,
                 compound: Compound,
                 criteria: Criteria | Iterable[Criteria] = None,
                 all_criteria: bool = False,
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
        super().__init__(compound.retention_time, None, criteria, all_criteria)
        self.compound = compound

    def create_peak(self, peak: PeakBounded):
        return PeakCompound(peak.parent, peak.bounds, self.compound, peak.id_)


class PeakCompound(PeakBounded):
    def __init__(self, parent: PeakParent, bounds: slice, compound: Compound, id_: int = None):
        super().__init__(parent, bounds, id_)
        self.compound = compound

    def get_stats(self) -> OrderedDict:
        dict_ = self.compound.get_stats()
        dict_.update(super().get_stats())
        return dict_
