from dataclasses import dataclass
import pathlib
from datetime import datetime, timedelta

from chem_analysis.nmr.parameters import NMRParameters


@dataclass(slots=True)
class NMRParameters(NMRParameters):
    is2D: bool = None
    sample_name: str = None
    solvent: str = None
    date_start: datetime = None
    nucleus: str = None
    temperature_coil: float = None  # Kelvin
    software_version: str = None
    software: str = None
    instrument: str = None
    pulse_sequence: str = None
    pulse_angle: int = None
    pulse_length: timedelta = None
    number_scans: int = None
    number_points: int = None
    spectrometer_frequency: float = None
    repetition_time: timedelta = None
    dwell_time: timedelta = None
    receiver_gain: float = None

    # acquisition_time: timedelta = None
    sizeTD1: int = None
    spectral_width: float = None
    sizeTD2: int = None
    carrier: float = None
    probe: str = None
    instrument_position: int = None
    shift_points: int = None

    # @property
    # def dwell_time(self) -> float:
    #     return 1 / self.spectral_width

    @property
    def acquisition_time(self) -> timedelta:
        return self.repetition_time - self.dwell_time


parse_mapping = {
    "Solvent": "solvent",
    "Sample": "sample_name",
    "Custom": None,
    "startTime": {"label": "date_start", "func": lambda x: datetime.fromisoformat(x)},
    "PreProtocolDelay": None,
    "acqDelay": None,  # shift_points ??
    "b1Freq": "spectrometer_frequency",
    "bandwidth": None,
    "dwellTime": {"label": "dwell_time", "func": lambda x: timedelta(milliseconds=x)},
    "experiment": {"label": "is2D", "func": lambda x: x == "2D"},
    "expName": None,
    "nrPnts": "number_points",
    "nrScans": "number_scans",
    "repTime": {"label": "repetition_time", "func": lambda x: timedelta(milliseconds=x)},
    "rxChannel": "nucleus",
    "rxGain": "receiver_gain",
    "lowestFrequency": None,
    "totalAcquisitionTime": None,
    "graphTitle": None,
    "linearPrediction": None,
    "userData": None,
    "90Amplitude": None,
    "pulseLength": {"label": "pulse_length", "func": lambda x: timedelta(milliseconds=x)},
    "pulseAngle": "pulse_angle",
    "ComputerName": None,
    "UserName": None,
    "SpinsolveUser": None,
    "ProtocolDataID": None,
    "Protocol": "pulse_sequence",
    "Options": None,
    "Spectrometer": None,
    "InstrumentType": "instrument",
    "InstrumentCode": None,
    "Software": "software_version",
    "WindowsLoggedUser": "software",
    "BackupLocation": None,
    "Shim_Timestamp": None,
    "Shim_Width50": None,
    "Shim_Width055": None,
    "Shim_SNR": None,
    "Shim_ReferencePeakIndex": None,
    "Shim_ReferencePPM": None,
    "Reference_Timestamp": None,
    "Reference_File": None,
    "StartProtocolTemperatureMagnet": None,
    "StartProtocolTemperatureBox": None,
    "StartProtocolTemperatureRoom": None,
    "CurrentTemperatureMagnet": {"label": "temperature_coil", "func": lambda x: x+273.15},
    "CurrentTemperatureBox": None,
    "CurrentTemperatureRoom": None,
    "CurrentTime": None,
}


def parse_acqu_file(path: pathlib.Path) -> NMRParameters:
    """

    Parameters
    ----------
    path:
        should finish with "/acqu.par"

    Returns
    -------

    """
    with open(path / "acqu.par", mode='r') as f:
        text = f.read()

    lines = text.strip().split("\n")
    parameters = NMRParameters()

    for line in lines:
        parameter, value = line.split("=")
        parameter = parameter.strip()
        value = value.strip()
        if '"' in value:
            value = value.replace('"', "")
        elif value[0].isdigit():
            value = float(value)
            if value == int(value):
                value = int(value)

        parameter = parse_mapping[parameter]
        if parameter is None:
            continue

        if isinstance(parameter, dict):
            value = parameter["func"](value)
            parameter = parameter["label"]

        setattr(parameters, parameter, value)

    return parameters


def main():
    import os
    MYDIR = r"C:\Users\nicep\Desktop\DW3-4"
    folder_paths = []
    for entry_name in os.listdir(MYDIR):
        entry_path = os.path.join(MYDIR, entry_name)
        if os.path.isdir(entry_path):
            folder_paths.append(entry_path)

    folder_paths = sorted(folder_paths, key=lambda x: int(x.split("_")[1]))
    for folder in folder_paths:
        folder_ = pathlib.Path(folder)
        result = parse_acqu_file(folder_)
        print(int(folder.split("_")[1]), result.date_start)


if __name__ == "__main__":
    main()
