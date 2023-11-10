import enum
from dataclasses import dataclass
from datetime import datetime

import numpy as np


class NMRExperiments(str, enum.Enum):
    PROTON = "1H"
    TOCSY = "TOCSY"
    NOESY = "NOESY"
    COSY = "COSY"
    CARBON = "13C"
    DEPT = "DEPT"
    HETCOR = "HETCOR"
    HSQC = "HSQC"
    HSBC = "HSBC"
    FLUORINE = "19F"
    PHOSPHORUS = "31P"


@dataclass(slots=True)
class NMRParameters:
    """

    """
    type_: NMRExperiments = None
    solvent: str = None
    date_start: datetime = None
    nmr_software_version: str = None
    spectrometer_frequency: float = None
    temperature: float = None  # Kelvin

    number_scans: int = None
    repetition_time: float = None
    pulse_angle: float = None  # degree
    pulse_length: float = None
    acquisition_time: float = None
    number_points: float = None


    sizeTD1: int = None
    spectral_width: float = None
    sizeTD2: int = None


    instrument: str = None
    pulse_sequence: str = None
    probe: str = None
    receiver_gain: float = None

    spectral_size: int = None

    instrument_position: int = None
    shift_points: int = None

    def compute_time(self):
        return np.linspace(
                0,
                (self.parameters.sizeTD2 - 1) * self.parameters.dwell_time,
                self.parameters.sizeTD2
            )

    @property
    def dwell_time(self) -> float:
        return 1 / self.spectral_width


class NMRParameters2D(NMRParameters):
    pass



"""
Acquisition Time (AT): The acquisition time, often denoted as AT, is the total duration for which the FID signal is recorded. It determines the length of the time domain signal. AT is usually set by the user and depends on the desired spectral width and resolution of the NMR spectrum.

Number of Data Points (N): The number of data points or data values collected during the FID acquisition is determined by the user and should be specified before the experiment. It affects the digital resolution and sensitivity of the resulting NMR spectrum. The greater the number of data points, the finer the spectral resolution, but the longer the acquisition time.

Dwell Time (DT): The dwell time, often denoted as DT, is the time interval between each data point collected in the FID. It is the reciprocal of the spectral width (SW) and is calculated as DT = 1/SW. DT determines the time resolution of the FID.

Zero Filling (ZF): Zero filling is a processing step performed after data acquisition to increase the number of data points in the FID. Zero filling involves interpolating zeros between the acquired data points to increase the apparent resolution. It can be used to enhance the resolution without changing the acquisition time.

Apodization Function: Apodization is a mathematical operation applied to the FID to modify its shape and reduce sideband artifacts in the Fourier-transformed spectrum. Common apodization functions include exponential multiplication, Gaussian multiplication, and sine-bell functions. The choice of apodization function affects the line shape and resolution of the resulting spectrum.

Frequency Offset (O1): The frequency offset (O1) represents the position of the spectral center relative to the transmitter frequency. It is set according to the region of interest in the NMR spectrum. The frequency offset can be applied during the acquisition to shift the FID to the desired spectral region.
"""