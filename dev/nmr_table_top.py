import glob
import os
import pathlib
import time

import numpy as np

from chem_analysis.nmr.parse_spinsolve_parameters import parse_acqu_file


def process_one(path: pathlib.Path) -> tuple[float, np.ndarray]:
    parameters = parse_acqu_file(path)

    data = np.loadtxt(path / "spectrum_processed.csv", delimiter=",", skiprows=1)

    return parameters.date_start.timestamp(), data


def process_many(path: pathlib.Path):
    files = tuple(os.scandir(path))
    files = sorted(files, key=lambda x: int(x.name[1:]))
    times = np.empty(len(files))
    ppm = None
    data = None
    ppm_slice = slice(-2, 12)  # ppm
    ppm_slice_index = None
    for i, file in enumerate(files):
        t, d = process_one(path / file.name)

        if ppm is None:
            ppm = d[:, 0]
            ppm_slice_index = slice(np.argmin(np.abs(ppm-ppm_slice.start)), np.argmin(np.abs(ppm-ppm_slice.stop)))
            ppm = ppm[ppm_slice_index]
            data = np.empty((len(files), len(ppm)), dtype=d.dtype)
        data[i, :] = d[ppm_slice_index, 1]
        times[i] = t

    return times, ppm, data


def repackage(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    data = np.empty((len(y)+1, len(x)+1), dtype=z.dtype)
    data[0, 1:] = x
    data[1:, 0] = y
    data[1:, 1:] = z
    return data


import pyarrow as pa

def write(arr, name):
    arrays = [pa.array(col) for col in arr]
    names = [str(i) for i in range(len(arrays))]
    batch = pa.RecordBatch.from_arrays(arrays, names=names)
    with pa.OSFile(name, 'wb') as sink:
        with pa.RecordBatchStreamWriter(sink, batch.schema) as writer:
            writer.write_batch(batch)


def read(name):
    source = pa.memory_map(name, 'r')
    table = pa.ipc.RecordBatchStreamReader(source).read_all()
    data = np.empty((table.num_columns, table.num_rows))
    for col in range(table.num_columns):
        data[col, :] = table.column(str(col))
    return data


def rename_spectra(path: pathlib.Path):
    folders = tuple(os.scandir(path))
    for folder in folders:
        modified_lines = []
        with open(pathlib.Path(folder.path) / "acqu.par", 'r') as file:
            for line in file:
                if line.startswith("Sample"):
                    # Replace the line that starts with "Sample ="
                    modified_lines.append(line.replace("DW2-5-1", f"DW2-5-1_{folder.name[1:]}"))
                else:
                    modified_lines.append(line)

        # Open the file in write mode and write the modified lines
        with open(pathlib.Path(folder.path) / "acqu.par", 'w') as file:
            file.writelines(modified_lines)


def main():
    path = pathlib.Path(r"C:\Users\nicep\Desktop\DW2-7")
    times, ppm, data = process_many(path)
    data = repackage(ppm, times, data)

    # csv
    # start = time.perf_counter()
    # np.savetxt("DW2_5_1_NMR.csv", data, delimiter=",")
    # print(time.perf_counter() - start)
    #
    # start = time.perf_counter()
    # data = np.loadtxt("DW2_5_1_NMR.csv", delimiter=",")
    # print(time.perf_counter() - start)

    # feather
    start = time.perf_counter()
    write(data, "DW2_7_NMR.feather")
    print(time.perf_counter()-start)
    #
    # start = time.perf_counter()
    # data_ = read("DW2_5_1.feather")
    # print(time.perf_counter()-start)
    # print(data_.shape)


def main2():
    path = pathlib.Path(r"C:\Users\nicep\Desktop\DW2-5-1")
    rename_spectra(path)


if __name__ == "__main__":
    main()
    # main2()
