
from chem_analysis.utils.printing_tables import StatsTable
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D


class ResultPeakPicking:
    def __init__(self, signal: Signal):
        self.signal = signal
        self.indexes = None

    def __str__(self):
        return f"# of Peaks: {len(self.indexes)}"

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        return iter(self.indexes)

    def __len__(self):
        return len(self.indexes)

    def values(self):
        return self.signal.x[self.indexes]

    def stats(self) -> StatsTable:
        headers = ["peak #", "index", "x", "y"]
        rows = []
        for i, index_ in enumerate(self.indexes):
            rows.append([i, index_, self.signal.x[index_], self.signal.y[index_]])
        return StatsTable(rows, headers)

    def add_peak(self, peak: int):
        self.indexes.append(peak)


class ResultPeakPickingArray:
    def __init__(self, signal2d: Signal2D):
        self.signal2d = signal2d
        self.results: list[ResultPeakPicking] = []

    def __iter__(self):
        return iter(self.results)

    def __len__(self):
        return len(self.results)

    def add_result(self, result: ResultPeakPicking):
        self.results.append(result)

    def stats(self) -> StatsTable:
        table = None
        for result in self.results:
            if table is None:
                table = result.stats()
            else:
                table.join(result.stats())

        return table
