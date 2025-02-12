from typing import Iterator
from collections import OrderedDict

from chem_analysis.analysis.peak import Peak
from chem_analysis.utils.printing_tables import StatsTable


class ResultPeaks:
    def __init__(self, peaks: list[Peak] | None = None) -> None:
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

    @property
    def name(self) -> str | None:
        if self.peaks:
            return str(self.peaks[0].parent.name)
        return None

    def to_dict(self) -> list[OrderedDict]:
        dicts_ = []
        for peak in self.peaks:
            dicts_.append(peak.properties.to_dict())

        return dicts_

    def stats_table(self) -> StatsTable:
        return StatsTable.from_list_dicts(self.to_dict())

    def add_peak(self, peak: Peak):
        self.peaks.append(peak)


class ResultPeaks2D:
    def __init__(self, results: list[ResultPeaks] = None):
        self.results = results or []

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
