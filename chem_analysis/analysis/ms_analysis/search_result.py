from typing import Iterator, OrderedDict

from chem_analysis.analysis.peak import PeakParent
from chem_analysis.analysis.integration.sec_peak import PeakIntegration
from chem_analysis.analysis.peak_result import ResultPeaks
from chem_analysis.gc_lc.library.compound import Compound


class PeakCompound(PeakIntegration):
    def __init__(self, parent: PeakParent, bounds: slice, compound: Compound, id_: int = None):
        super().__init__(parent, bounds, id_)
        self.compound = compound


class ResultCompoundSearch(ResultPeaks):
    def __init__(self, peaks: list[PeakCompound] = None):
        self.peaks = peaks or []

    def __str__(self):
        return f"# of Peaks: {len(self)}"

    def __repr__(self):
        return self.__str__()

    def __iter__(self) -> Iterator[Peak]:
        return iter(self.peaks)

    def __len__(self):
        return len(self.peaks)

    def __getitem__(self, item: int | slice):
        return self.peaks[item]

    def get_stats(self) -> list[OrderedDict]:
        dicts_ = []
        for peak in self.peaks:
            dicts_.append(peak.get_stats())

        return dicts_

    def add_peak(self, peak: Peak):
        self.peaks.append(peak)