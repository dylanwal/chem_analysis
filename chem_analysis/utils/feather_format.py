import pathlib
from typing import Sequence

import numpy as np
import pyarrow as pa


def numpy_to_feather(array_: np.ndarray, file_path: str | pathlib.Path):
    """
    Save numpy array to feather file

    Parameters
    ----------
    array_:
        numpy array
    file_path:
        file path

    """
    if not isinstance(file_path, pathlib.Path):
        file_path = pathlib.Path(file_path)
    if file_path.suffix != ".feather":
        file_path = file_path.with_suffix(".feather")

    arrays = [pa.array(col) for col in array_]
    names = [str(i) for i in range(len(arrays))]
    batch = pa.RecordBatch.from_arrays(arrays, names=names)
    with pa.OSFile(str(file_path), 'wb') as sink:
        with pa.RecordBatchStreamWriter(sink, batch.schema) as writer:
            writer.write_batch(batch)


def feather_to_numpy(file_path: str | pathlib.Path) -> np.ndarray:
    """

    Parameters
    ----------
    file_path:
        feather file

    Returns
    -------

    """
    if isinstance(file_path, pathlib.Path):
        file_path = str(file_path)

    source = pa.memory_map(file_path, 'r')
    table = pa.ipc.RecordBatchStreamReader(source).read_all()
    data = np.empty((table.num_columns, table.num_rows))
    for col in range(table.num_columns):
        data[col, :] = table.column(str(col))
    return data


def pack_time_series(x: np.ndarray, time_: np.ndarray, z: np.array) -> np.ndarray:
    data = np.empty((len(time_) + 1, len(x) + 1), dtype=z.dtype)
    data[0, 0] = 0
    data[0, 1:] = x
    data[1:, 0] = time_
    data[1:, 1:] = z
    return data


def unpack_time_series(data: np.array) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[0, 1:]
    time_ = data[1:, 0]
    z = data[1:, 1:]
    return x, time_, z


def unpack_and_merge_time_series_feather_files(paths: Sequence[str | pathlib.Path]) \
        -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if len(paths) == 0:
        raise ValueError("No files to unpack.")

    z_data = []
    time_data = []
    for path in paths:
        x, time_, z = unpack_time_series(feather_to_numpy(path))
        z_data.append(z)
        time_data.append(time_)

    return x, np.concatenate(time_data), np.concatenate(z_data)
