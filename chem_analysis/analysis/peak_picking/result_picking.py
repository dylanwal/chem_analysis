
from chem_analysis.analysis.peak import PeakDiscrete
from chem_analysis.analysis.peak_result import ResultPeaks, ResultPeaks2D
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D


class PeakPicking(PeakDiscrete):
    def __init__(self, parent: Signal, index: int, id_: int = None):
        super().__init__(parent, index, id_)

    @property
    def signal(self) -> Signal:
        return self.parent


class ResultPicking(ResultPeaks):
    def __init__(self, peaks: list[PeakPicking] = None, signal: Signal = None) -> None:
        super().__init__(peaks)
        self.signal = signal


class ResultPicking2D(ResultPeaks2D):
    def __init__(self, results: list[ResultPeaks] = None, signal: Signal2D = None) -> None:
        super().__init__(results)
        self.signal = signal
