import pathlib
import struct
import os
import configparser
from datetime import datetime

from pydantic import BaseModel

from chem_analysis.mass_spec.ms_parameters import MSParameters


class MSScanParameters(BaseModel):
    scan_start: float
    low_mass: int
    high_mass: int


class Timestamps(BaseModel):
    start: datetime
    acq_start: datetime
    perform_inj: datetime
    start_inj: datetime
    start_run: datetime
    start_post_run: datetime
    start_post_acq: datetime
    meth_done: datetime


class MassHunterParameters(MSParameters):
    def __init__(self,
                 MS_scan: MSScanParameters = None,
                 method: pathlib.Path = None,
                 date_times: Timestamps = None,
                 ):
        super().__init__()
        self.MS_scan = MS_scan
        self.method_path = method
        self.date_times = date_times


def parse_D_folder(folder_path: str):
    parameters = MassHunterParameters()

    # pre_post.ini
    file_path = os.path.join(folder_path, 'pre_post.ini')
    if os.path.exists(file_path):
        parse_pre_post_ini(parameters, file_path)

    # data.ms
    file_path = os.path.join(folder_path, 'data.ms')
    if os.path.exists(file_path):
        parse_data_ms(parameters, file_path)

    # FID1A.ch
    file_path = os.path.join(folder_path, 'FID1A.ch')
    if os.path.exists(file_path):
        parse_fid1a_ch(parameters, file_path)


def parse_pre_post_ini(parameters: MassHunterParameters, file_path: str):
    config = configparser.ConfigParser()
    config.read(file_path)

    if 'Scan Parameters' in config:
        parameters.MS_scan = MSScanParameters(
            scan_start=config['Scan Parameters']['scanstart1'],
            low_mass=config['Scan Parameters']['lowmass1'],
            high_mass=config['Scan Parameters']['highmass1']
        )

    if 'Sequence' in config:
        parameters.method_path = pathlib.Path(config['Sequence']['_methpath$'] + config['Sequence']['_methfile$'])

    if 'Timings' in config:
        parameters.date_times = Timestamps(
            start=config['Timings']['TS1_______Start$'],
            acq_start=config['Timings']['TS2____AcqStart$'],
            perform_inj=config['Timings']['TS3__PerformInj$'],
            start_inj=config['Timings']['TS4____StartInj$'],
            start_run=config['Timings']['TS5____StartRun$'],
            start_post_run=config['Timings']['TS6StartPostRun$'],
            start_post_acq=config['Timings']['TS7StartPostAcq$'],
            meth_done=config['Timings']['TS8____MethDone$'],
        )


def parse_data_ms(parameters: MSParameters, file_path: str):
    with open(file_path, 'rb') as f:
        while True:
            try:
                # Read text (assuming it's a null-terminated string)
                text = b''
                while True:
                    char = f.read(1)
                    text += char
                text = text.decode('utf-8')

                # Read integer
                integer = struct.unpack('i', f.read(4))[0]

                print(f'Text: {text}, Integer: {integer}')
            except struct.error:
                break  # Reached end of file

        # head = read_string(f, offset=0, gap=1)
        # if head in ['179', '181']:
        #     return parse_ch_fid(path, head)
        # elif head in ['130', '30']:
        #     return parse_ch_other(path, head)


def parse_fid1a_ch(file_path: str):
    pass
    # with open(path, 'rb') as f:
    #     head = read_string(f, offset=0, gap=1)
    #     if head in ['179', '181']:
    #         return parse_ch_fid(path, head)
    #     elif head in ['130', '30']:
    #         return parse_ch_other(path, head)
    #     return None


def read_string(f, offset, gap=2):
    """
    Extracts a string from the specified offset.

    This method is primarily useful for retrieving metadata.

    Args:
        f (_io.BufferedReader): File opened in 'rb' mode.
        offset (int): Offset to begin reading from.
        gap (int): Distance between two adjacent characters.

    Returns:
        String at the specified offset in the file header.

    """
    f.seek(offset)
    str_len = struct.unpack("<B", f.read(1))[0] * gap
    try:
        return f.read(str_len)[::gap].decode().strip()
    except Exception:
        return ""


def main():
    path = r"C:\Users\nicep\Desktop\research_wis\data\10\11\10_11\DJW-10-11-2h-batch-PPh3.D"
    result = parse_D_folder(path)
    print(result)


if __name__ == "__main__":
    main()
