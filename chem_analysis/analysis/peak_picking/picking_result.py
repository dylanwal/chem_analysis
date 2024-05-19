from typing import Iterator
from collections import OrderedDict

from chem_analysis.analysis.peak import Peak
from chem_analysis.utils.printing_tables import StatsTable
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D


class ResultPeaks:
    def __init__(self, signal: Signal):
        self.signal = signal
        self.peaks: list[Peak] = []

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

    # def values(self):
    #     return self.signal.x[self.peaks]

    def get_stats(self) -> list[OrderedDict]:
        dicts_ = []
        for peak in self.peaks:
            dicts_.append(peak.get_stats())

        return dicts_

    def stats_table(self) -> StatsTable:
        return StatsTable.from_list_dicts(self.get_stats())

    def add_peak(self, peak: Peak):
        self.peaks.append(peak)


class ResultPeakArray:
    def __init__(self, signal2d: Signal2D):
        self.signal2d = signal2d
        self.results: list[ResultPeaks] = []

    def __iter__(self):
        return iter(self.results)

    def __getitem__(self, item: int | slice):
        return self.results[item]

    def __len__(self):
        return len(self.results)

    def add_result(self, result: ResultPeaks):
        self.results.append(result)

    def get_stats(self) -> StatsTable:
        table = None
        for result in self.results:
            if table is None:
                table = result.stats_table()
            else:
                table.join(result.stats_table())

        return table
