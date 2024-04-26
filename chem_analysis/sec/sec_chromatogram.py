
from chem_analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.sec.sec_signal import SECSignal
from chem_analysis.sec.sec_calibration import SECCalibration


class SECChromatogram(Chromatogram):

    def __init__(self, data:  list[SECSignal], calibration: SECCalibration = None):
        self.calibration = calibration

        super().__init__(data)
