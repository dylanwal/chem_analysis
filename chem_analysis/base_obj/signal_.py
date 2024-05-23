from typing import Sequence
import pathlib

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.processing.processor import Processor


def validate_input(x_raw: np.ndarray, data_raw: np.ndarray):
    if len(x_raw.shape) != 1:
        raise ValueError(f"'x_raw' must shape 1. \n\treceived: {x_raw.shape}")
    if len(data_raw.shape) != 1:
        raise ValueError(f"'data_raw' must shape 1. \n\treceived: {data_raw.shape}")
    if x_raw.shape != data_raw.shape:
        raise ValueError(f"'x_raw' and 'data_raw' must have same shape. \n\treceived: x_raw:{x_raw.shape} || data_raw: "
                         f"{data_raw.shape}")


class Signal:
    """ signal

    A signal is any x-y data.

    """
    __count = 0

    def __init__(self,
                 x_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        """

        Parameters
        ----------
        x_raw: np.ndarray[i]
            raw x data, i length
        data_raw: np.ndarray[i]
            raw y data, i length
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        name: str
            user defined name
        """
        validate_input(x_raw, data_raw)

        self.x_raw = x_raw
        self.data_raw = data_raw
        self.id_ = id_ or Signal.__count
        Signal.__count += 1
        self.name = name or f"signal_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"

        self.processor = Processor()
        self._x = None
        self._data = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.x)

    def _process(self):
        self._x, self._data = self.processor.run(self.x_raw, self.data_raw)

    @property
    def x(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._x

    @property
    def data(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._data

    @property
    def y(self) -> np.ndarray:
        return self.data

    def y_normalized_by_max(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return self.y/np.max(self.y)
        return general_math.normalize_by_max_with_x_range(x=self.x, y=self.y, x_range=x_range)

    def y_normalized_by_area(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return general_math.normalize_by_area(x=self.x, y=self.y)
        return general_math.y_normalized_by_area_with_x_range(x=self.x, y=self.y, x_range=x_range)

    def to_dict(self, sanitize: bool = False) -> dict:
        dict_ = {
            "name": self.name,
            "x_label": self.x_label,
            "y_label": self.y_label,
            "id_": self.id_
        }
        if sanitize:
            dict_["data"] = np.column_stack([self.x, self.y]).tolist()
        else:
            dict_["data"] = np.column_stack([self.x, self.y])

        return dict_

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################
    def to_json(self, path: str | pathlib.Path, encoding: str = "utf-8", **kwargs):
        import json

        kwargs = kwargs or dict()
        with open(path, 'w', encoding=encoding) as file:
            json.dump(self.to_dict(sanitize=True), file, **kwargs)

    def to_csv(self, path: str | pathlib.Path, headers: bool = False, encoding: str = "utf-8"):
        kwargs = {"encoding": encoding}
        if headers:
            kwargs["headers"] = [self.x_label, self.y_label]
        np.savetxt(path, np.column_stack((self.x, self.y)), delimiter=",", **kwargs)

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        np.save(path, np.column_stack((self.x, self.y)), **kwargs)

    def to_feather(self, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import numpy_to_feather

        headers = [self.x_label, self.y_label]
        numpy_to_feather(np.column_stack((self.x, self.y)), path, headers=headers)

    def to_parquet(self, path: str | pathlib.Path):
        import pyarrow as pa
        import pyarrow.parquet as pq
        # TODO:
        raise NotImplementedError()
        arrays = [pa.array(self.x), pa.array(self.y)]
        table = pa.Table().from_arrays(arrays=arrays, names=("x", "y"))
        pq.write_table(table, path)

    @classmethod
    def from_parquet(cls, path: str | pathlib.Path, **kwargs):
        if isinstance(path, str):
            path = pathlib.Path(path)

        # TODO:
        raise NotImplementedError()

    @classmethod
    def from_json(cls, path: str | pathlib.Path, encoding: str = "utf-8", **kwargs):
        if isinstance(path, str):
            path = pathlib.Path(path)

        import json
        with open(path, 'r', encoding=encoding) as file:
            data = json.load(file)

        signal = np.array(data['data'])
        x = signal[:, 0]
        y = signal[:, 1]
        return cls(x, y, x_label=data['x_label'], y_label=data['y_label'], id_=data['id_'])

    @classmethod
    def from_csv(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        x, y, x_label, y_label = load_csv(path)
        return cls(x, y, x_label=x_label, y_label=y_label)

    @classmethod
    def from_npy(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        x, y = np.load(str(path))
        return cls(x, y)

    @classmethod
    def from_feather(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        from chem_analysis.utils.feather_format import feather_to_numpy
        data, headers = feather_to_numpy(path)
        x, y = data[:, 0], data[:, 1]
        if headers[0] != "0":
            x_label = headers[0]
            y_label = headers[1]
        else:
            x_label = y_label = None

        return cls(x, y, x_label=x_label, y_label=y_label)


def load_csv(path: pathlib) -> tuple[np.ndarray, np.ndarray, str | None, str | None]:
    import csv

    data = []
    x_label = None
    y_label = None

    with open(path, 'r') as file:
        csv_reader = csv.reader(file)

        first_row = next(csv_reader, None)
        if len(first_row) != 2:
            raise ValueError("Data not correct format.")
        if any(cell.isalpha() for cell in first_row):
            # If the first row contains non-numeric values, consider it as column labels
            x_label, y_label = first_row
        else:
            # If the first row contains numbers, treat them as data and set labels to None
            data.append([float(cell) for cell in first_row])

        for row in csv_reader:
            data.append([float(cell) for cell in row])

    data_array = np.array(data)
    return data_array[:, 0], data_array[:, 1], x_label, y_label
