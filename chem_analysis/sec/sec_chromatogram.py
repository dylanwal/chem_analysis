
from chem_analysis.base_obj.calibration import Calibration
from chem_analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.sec.sec_signal import SECSignal


class SECChromatogram(Chromatogram):

    def __init__(self, data: Union[pd.DataFrame, Signal, list[Signal]], calibration: Calibration = None):
        if not isinstance(cal, Cal):
            cal = Cal(cal)
        self.calibration = calibration

        super().__init__(data)
