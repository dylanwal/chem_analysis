from __future__ import annotations
from typing import Iterable
from dataclasses import dataclass
import pathlib
from datetime import datetime
import math


@dataclass(slots=True)
class NMRParameters:
    is2D: bool = None
    sizeTD1: int = None
    spectral_width: float = None
    sizeTD2: int = None
    carrier: float = None
    date_start: datetime = None
    number_scans: int = None
    solvent: str = None
    temperature_coil: float = None  # Kelvin
    instrument: str = None
    pulse_sequence: str = None
    probe: str = None
    receiver_gain: float = None
    spectrometer_frequency: float = None
    spectral_size: int = None
    top_spin_version: str = None
    instrument_position: int = None
    shift_points: int = None
    nucleus: str = None

    @property
    def dwell_time(self) -> float:
        return 1 / self.spectral_width


parseDict = {3: "i4", 4: "f8"}


def get_line(lines: Iterable[str], startswith: str) -> str | None:
    for line in lines:
        if line.startswith(startswith):
            return line

    return None


def parse_acqus_file(path: pathlib.Path) -> NMRParameters:
    """

    Parameters
    ----------
    path:
        should finish with "/acqus"

    Returns
    -------

    """
    with open(path, mode='r') as f:
        text = f.read()

    lines = text.split("\n")
    acqusfile = NMRParameters()

    line = get_line(lines, "##$SW_h=").strip("##$SW_h= ")
    acqusfile.spectral_width = float(line)

    line = get_line(lines, "##$TD").strip("##$TD= ")
    acqusfile.sizeTD2 = int(int(line) / 2)

    line = get_line(lines, "##$TD").strip("##$TD= ")
    acqusfile.carrier = float(line) * 1e6

    line = get_line(lines, "##$DATE_START").strip("##$DATE_START= ")
    acqusfile.date_start = datetime.fromtimestamp(int(line))

    line = get_line(lines, "##$NS=").strip("##$NS= ")
    acqusfile.number_scans = int(line)

    line = get_line(lines, "##$SOLVENT=").strip("##$SOLVENT= ").replace("<", "").replace(">", "")
    acqusfile.solvent = line

    line = get_line(lines, "##$shimCoilTempK=").strip("##$shimCoilTempK= ")
    acqusfile.temperature_coil = float(line)

    line = get_line(lines, "##$INSTRUM=").strip("##$INSTRUM= ").replace("<", "").replace(">", "")
    acqusfile.instrument = line

    line = get_line(lines, "##$PULPROG=").strip("##$PULPROG= ").replace("<", "").replace(">", "")
    acqusfile.pulse_sequence = line

    line = get_line(lines, "##$PROBHD=").strip("##$PROBHD= ").replace("<", "").replace(">", "")
    acqusfile.probe = line

    line = get_line(lines, "##$RG=").strip("##$RG= ")
    acqusfile.receiver_gain = float(line)

    line = get_line(lines, "##$SFO1=").strip("##$SFO1= ")
    acqusfile.spectrometer_frequency = float(line)

    line = get_line(lines, "##$TD=").strip("##$TD= ")
    acqusfile.spectral_size = float(line)

    line = get_line(lines, "##TITLE").strip("##TITLE= Parameter file, TopSpin ")
    acqusfile.top_spin_version = line

    line = get_line(lines, "##$HOLDER").strip("##$HOLDER= ")
    acqusfile.instrument_position = int(line)

    line = get_line(lines, "##$GRPDLY").strip("##$GRPDLY= ")
    acqusfile.shift_points = int(math.floor(float(line)))

    line = get_line(lines, "##$NUC1").strip("##$NUC1= ").replace("<", "").replace(">", "")
    acqusfile.nucleus = line

    return acqusfile
