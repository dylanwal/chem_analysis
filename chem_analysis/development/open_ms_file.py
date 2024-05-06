import base64
import struct

import numpy as np

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


def parse_fid(file_path: str):
    data = {}
    with open(file_path, 'br') as file:
        data['version'] = fpascal(file, 0)

        data["file_info"] = fpascal(file, 4, 'UTF-8')
        data["sample_name"] = fpascal(file, 24, 'UTF-8')
        data["sample_info"] = fpascal(file, 86, 'UTF-8')
        data["operator"] = fpascal(file, 148, 'UTF-8')
        data["date_time"] = fpascal(file, 178, 'UTF-8')
        data["instmodel"] = fpascal(file, 208, 'UTF-8')
        data["inlet"] = fpascal(file, 218, 'UTF-8')
        data["method_name"] = fpascal(file, 228, 'UTF-8')
        data["seqindex"] = fnumeric(file, 252, '>h')
        data["vial"] = fnumeric(file, 254, '>h')
        data["replicate"] = fnumeric(file, 256, '>h')
        data["start_time"] = fnumeric(file, 282, '>i') / 6E4  # min
        data["end_time"] = fnumeric(file, 286, '>i') / 6E4  # min
        data["channel_max"] = fnumeric(file, 290, '>i')
        data["channel_min"] = fnumeric(file, 294, '>i')
        data["dir_type"] = fnumeric(file, 258, '>h')
        data["dir_offset"] = fnumeric(file, 260, '>i') * 2 - 2
        data["data_offset"] = fnumeric(file, 264, '>i') * 2 - 2
        data["num_records"] = fnumeric(file, 278, '>i')
        data['glp_flag'] = fnumeric(file, 3085, '>i')
        data["data_source"] = fpascal(file, 3089, 'UTF-16')
        data["firmware_rev"] = fpascal(file, 3601, 'UTF-16')
        data["software_rev"] = fpascal(file, 3802, 'UTF-16')
        data["intensity_units"] = "counts"
        data["channel_units"] = "m/z"
        data["time_units"] = "min"

        spectrum_offset, retention_time, total_abundance = fdirectory(file, data["dir_offset"], data['num_records'])
        data['intensity'] = total_abundance
        data['time'] = retention_time
        data['data_offset'] = spectrum_offset

        data['data'] = fscan(file, data['data_offset'])

    return data


def fpascal(f, offset, encoding: str = 'UTF-8') -> str:
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


def fnumeric(f, offset: int, encoding: str = '<B') -> int | float:
    f.seek(offset, 0)
    return struct.unpack(encoding, f.read(DATA_TYPE_SIZE[encoding[1]]))[0]


def fscan(f, offset):
    n = []
    mz = []
    intensity = []

    for i in range(len(offset)):
        f.seek(offset[i]+12, 0)
        n.append(struct.unpack('>h', f.read(2))[0])
        f.seek(4, 1)
        d = np.fromfile(f, dtype=np.dtype('uint16').newbyteorder('>'), count=2*n[i])
        mz.append(np.round(np.flip(d[::2]) / 20).astype('uint16'))
        intensity.append(np.flip(d[1::2]))

    # pack data
    data = np.zeros((len(mz), 1000), dtype="int16")
    for i, (m, y) in enumerate(zip(mz, intensity)):
        data[i, m] = y

    return data


def fdirectory(f, offset: int, num_records: int):
    f.seek(offset, 0)

    # # Read directory contents
    data = np.fromfile(f, dtype=np.dtype('int32').newbyteorder('>'), count=num_records*3)
    spectrum_offset = data[::3]
    retention_time = data[1::3]
    total_abundance = data[2::3]

    # Apply correction factors
    spectrum_offset = spectrum_offset * 2 - 2
    retention_time = retention_time / 6E4

    return spectrum_offset, retention_time, total_abundance


def main():
    file_path = r"C:\Users\nicep\Desktop\research_wis\data\10\10_13\DJW-10-13-600min-TMS.D\data.ms"
    data = parse_fid(file_path)

    x = data["time"]
    y = np.sum(data['data'], axis=1)
    import plotly.graph_objs as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y))
    fig.show()
    print(data)



if __name__ == "__main__":
    main()