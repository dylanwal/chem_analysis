import os
import configparser
import pathlib
import struct
from datetime import datetime
from typing import Any, BinaryIO
import logging

import numpy as np

logger = logging.getLogger(__name__)

# https://docs.python.org/3/library/struct.html#format-characters
DATA_TYPE_SIZE = {
    "c": 1,
    "b": 1,
    "B": 1,
    "?": 1,
    "h": 2,
    "H": 2,
    "i": 4,
    "I": 4,
    "l": 4,
    "L": 4,
    "q": 8,
    "Q": 8,
    "e": 2,
    "f": 4,
    "d": 8
}


def parse_D_folder(folder_path: pathlib.Path) -> tuple[dict | None, dict | None, dict | None]:
    ini_dict, ms_dict, fid_dict = None, None, None

    ini_path = folder_path / 'pre_post.ini'
    if ini_path.exists():
        ini_dict = parse_pre_post_ini(ini_path)

    ms_path = folder_path / 'data.ms'
    if ms_path.exists():
        ms_dict = parse_gcms(ms_path)

    fid_path = folder_path / 'FID1A.ch'
    if fid_path.exists():
        fid_dict = parse_fid(fid_path)

    return ini_dict, ms_dict, fid_dict


def parse_pre_post_ini(file_path: str | pathlib.Path) -> dict[str, Any]:
    config = configparser.ConfigParser()
    config.read(file_path)

    data = {}

    if 'Scan Parameters' in config:
        data["scan_start"] = float(config['Scan Parameters']['scanstart1'])
        data['low_mass'] = int(config['Scan Parameters']['lowmass1'])
        data['high_mass'] = int(config['Scan Parameters']['highmass1'])

    if 'Sequence' in config:
        data['method_path'] = pathlib.Path(config['Sequence']['_methpath$'] + config['Sequence']['_methfile$'])

    if 'Timings' in config:
        data['time_start'] = datetime.fromisoformat(config['Timings']['TS1_______Start$'])
        data['time_acq_start'] = datetime.fromisoformat(config['Timings']['TS2____AcqStart$'])
        data['time_perform_inj'] = datetime.fromisoformat(config['Timings']['TS3__PerformInj$'])
        data['time_start_inj'] = datetime.fromisoformat(config['Timings']['TS4____StartInj$'])
        data['time_start_run'] = datetime.fromisoformat(config['Timings']['TS5____StartRun$'])
        data['start_post_run'] = datetime.fromisoformat(config['Timings']['TS6StartPostRun$'])
        data['start_post_acq'] = datetime.fromisoformat(config['Timings']['TS7StartPostAcq$'])
        data['meth_done'] = datetime.fromisoformat(config['Timings']['TS8____MethDone$'])

    return data


def parse_fid(file_path: str | pathlib.Path):
    data = {}
    with open(file_path, 'br') as file:
        data['version'] = f_pascal(file, 0)

        data["file_info"] = f_pascal(file, 347, 'UTF-16')
        data["sample_name"] = f_pascal(file, 858, 'UTF-16')
        data["barcode"] = f_pascal(file, 1369, 'UTF-16')
        data["operator"] = f_pascal(file, 1880, 'UTF-16')
        data["date_time"] = datetime.strptime(f_pascal(file, 2391, 'UTF-16'), '%d %b %y  %I:%M %p')
        data["instmodel"] = f_pascal(file, 2492, 'UTF-16')
        data["inlet"] = f_pascal(file, 2533, 'UTF-16')
        data["method_name"] = f_pascal(file, 2574, 'UTF-16')
        data["seqindex"] = f_numeric(file, 252, '>h')
        data["vial"] = f_numeric(file, 254, '>h')
        data["replicate"] = f_numeric(file, 256, '>h')
        data["start_time"] = f_numeric(file, 282, '>f') / 6E4  # min
        data["end_time"] = f_numeric(file, 286, '>f') / 6E4  # min
        data["channel_max"] = f_numeric(file, 290, '>f')
        data["channel_min"] = f_numeric(file, 294, '>f')
        data["dir_type"] = f_numeric(file, 258, '>h')
        data["dir_offset"] = (f_numeric(file, 260, '>i') - 1) * 512
        data["data_offset"] = (f_numeric(file, 264, '>i') - 1) * 512
        data["num_records"] = f_numeric(file, 278, '>i')
        data['glp_flag'] = f_numeric(file, 3085, '>i')
        data["data_source"] = f_pascal(file, 3089, 'UTF-16')
        data["firmware_rev"] = f_pascal(file, 3601, 'UTF-16')
        data["software_rev"] = f_pascal(file, 3802, 'UTF-16')
        data["channel_detector"] = f_numeric(file, 4106, '>h')
        data["channel_units"] = f_pascal(file, 4172, 'UTF-16')
        data['channel_desc'] = f_pascal(file, 4213, 'UTF-16')

        # if f_numeric(file, data["data_offset"], '>h') == 2048:
        data["data_offset"] = data["data_offset"] + 2048

        data['data'] = f_double_array(file, data["data_offset"])
        data['time'] = np.linspace(data['start_time'], data['end_time'], data['data'].size)

    return data


def parse_gcms(file_path: str | pathlib.Path):
    data = {}
    with open(file_path, 'br') as file:
        data['version'] = int(f_pascal(file, 0))

        data["file_info"] = f_pascal(file, 4, 'UTF-8')
        data["sample_name"] = f_pascal(file, 24, 'UTF-8')
        data["sample_info"] = f_pascal(file, 86, 'UTF-8')
        data["operator"] = f_pascal(file, 148, 'UTF-8')
        data["date_time"] = datetime.strptime(f_pascal(file, 178, 'UTF-8'), '%d %b %y  %I:%M %p')
        data["instmodel"] = f_pascal(file, 208, 'UTF-8')
        data["inlet"] = f_pascal(file, 218, 'UTF-8')
        data["method_name"] = f_pascal(file, 228, 'UTF-8')
        data["seqindex"] = f_numeric(file, 252, '>h')
        data["vial"] = f_numeric(file, 254, '>h')
        data["replicate"] = f_numeric(file, 256, '>h')
        data["start_time"] = f_numeric(file, 282, '>i') / 6E4  # min
        data["end_time"] = f_numeric(file, 286, '>i') / 6E4  # min
        data["channel_max"] = f_numeric(file, 290, '>i')
        data["channel_min"] = f_numeric(file, 294, '>i')
        data["dir_type"] = f_numeric(file, 258, '>h')
        data["dir_offset"] = f_numeric(file, 260, '>i') * 2 - 2
        # data["data_offset"] = f_numeric(file, 264, '>i') * 2 - 2  # extracted below
        data["num_records"] = f_numeric(file, 278, '>i')
        data['glp_flag'] = f_numeric(file, 3085, '>i')
        data["data_source"] = f_pascal(file, 3089, 'UTF-16')
        data["firmware_rev"] = f_pascal(file, 3601, 'UTF-16')
        data["software_rev"] = f_pascal(file, 3802, 'UTF-16')
        data["intensity_units"] = "counts"
        data["channel_units"] = "m/z"
        data["time_units"] = "min"

        data_offset, time, intensity = f_directory(file, data["dir_offset"], data['num_records'])
        data['data_offset'], data['time'], data['intensity'] = data_offset, time, intensity
        data['data'] = f_scan(file, data['data_offset'])

    return data


def f_pascal(f, offset, encoding: str = 'UTF-8') -> str:
    f.seek(offset, 0)
    str_len = struct.unpack('<B', f.read(1))[0]
    if encoding == 'UTF-16':
        str_len *= 2
    bytes_ = f.read(str_len)
    str_ = bytes_.decode(encoding)

    if len(str_) > 512:
        str_ = ''
    else:
        str_ = str_.strip()

    return str_


def f_numeric(f: BinaryIO, offset: int, encoding: str = '<B') -> int | float:
    f.seek(offset, 0)
    return struct.unpack(encoding, f.read(DATA_TYPE_SIZE[encoding[1]]))[0]


def f_double_array(f: BinaryIO, offset: int) -> np.ndarray:
    f.seek(0, os.SEEK_END)
    file_size = f.tell()
    f.seek(offset, os.SEEK_SET)
    n = (file_size - offset) // 8
    data = np.fromfile(f, dtype=np.dtype("float64").newbyteorder('<'), count=n)

    return data


def f_scan(f: BinaryIO, offset: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    n = []
    mz = []
    intensity = []

    for i in range(len(offset)):
        f.seek(offset[i] + 12, 0)
        n.append(struct.unpack('>h', f.read(2))[0])
        f.seek(4, 1)
        d = np.fromfile(f, dtype=np.dtype('uint16').newbyteorder('>'), count=2 * n[i])
        mz.append(np.flip(d[::2]) / 20)
        intensity.append(np.flip(d[1::2]))

    # convert to int16 to int32 (some of int16 bits represent bit shifts)
    intensity_new = []
    for i, intensity_ in enumerate(intensity):
        e = np.bitwise_and(intensity_.astype('int32'), np.int32(49152))
        y = np.bitwise_and(intensity_.astype('int32'), np.int32(16383))

        while np.any(e != 0):
            y[e != 0] = np.left_shift(y[e != 0], 3)
            e[e != 0] = e[e != 0] - 16384

        intensity_new.append(y)

    return list(zip(mz, intensity_new))


def f_directory(f: BinaryIO, offset: int, num_records: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f.seek(offset, 0)

    # # Read directory contents
    data = np.fromfile(f, dtype=np.dtype('int32').newbyteorder('>'), count=num_records * 3)
    spectrum_offset = data[::3]
    retention_time = data[1::3]
    total_abundance = data[2::3]

    # Apply correction factors
    spectrum_offset = spectrum_offset * 2 - 2
    retention_time = retention_time / 6E4

    return spectrum_offset, retention_time, total_abundance
