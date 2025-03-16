from typing import Sequence
import pathlib
import copy

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.base_obj.parameters import Parameters


def check_array_inf_nan(x: np.ndarray, name: str):
    if np.any(np.isinf(x)):
        raise ValueError(f"The array '{name}' contains 'inf' values.")

    if np.any(np.isnan(x)):
        raise ValueError(f"The array '{name}' contains 'nan' values.")


def validate_input(x: np.ndarray, y: np.ndarray):
    if len(x.shape) != 1:
        raise ValueError(f"'x' must shape 1. \n\treceived: {x.shape}")
    if len(y.shape) != 1:
        raise ValueError(f"'y' must shape 1. \n\treceived: {y.shape}")
    if x.shape != y.shape:
        raise ValueError(f"'x' and 'y' must have same shape. \n\treceived: x:{x.shape} || y:{y.shape}")
    check_array_inf_nan(x, name="x")
    check_array_inf_nan(y, name="y")


class Signal:
    """ signal

    A signal is any x-y data.

    """
    __slots__ = "x", "y", "id_", "name", "x_label", "y_label", "parameters", "process_history", "extract_value"
    __count = 0

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 id_: int = None,
                 parameters: Parameters = None,
                 process_history: str | list[str] = False,
                 ):
        """

        Parameters
        ----------
        x: np.ndarray[i]
            raw x data, i length
        y: np.ndarray[i]
            raw y data, i length
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        name: str
            user defined name
        parameters: Parameters
            various meta-data
        process_history: list[str]
            List of processing methods the signal has been through
        """
        validate_input(x, y)
        x, y = general_math.check_for_flip(x, y)

        self.x = x
        self.y = y
        self.id_ = id_ or Signal.__count
        Signal.__count += 1
        self.name = name or f"signal_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"

        self.parameters = parameters
        self.extract_value = None

        if isinstance(process_history, str):
            process_history = [process_history]
        self.process_history = process_history or []

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self.x)})"
        return text

    def copy_with(self, x: np.ndarray, y: np.ndarray, deep: bool = True):
        if deep:
            copy_method = copy.deepcopy
        else:
            copy_method = copy.copy

        return Signal(
                x,
                y,
                name=copy_method(self.name),
                x_label=copy_method(self.x_label),
                y_label=copy_method(self.y_label),
                parameters=copy_method(self.parameters),
                process_history=copy_method(self.process_history)
            )

    def y_normalized_by_max(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return self.y / np.max(self.y)
        return general_math.normalize_by_max_with_x_range(x=self.x, y=self.y, x_range=x_range)

    def y_normalized_by_area(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return general_math.normalize_by_area(x=self.x, y=self.y)
        return general_math.y_normalized_by_area_with_x_range(x=self.x, y=self.y, x_range=x_range)

    def to_dict(self, data_as_list: bool = False) -> dict:
        dict_ = {
            "name": self.name,
            "x_label": self.x_label,
            "y_label": self.y_label,
            "id_": self.id_
        }
        if data_as_list:
            dict_["y"] = np.column_stack([self.x, self.y]).tolist()
        else:
            dict_["y"] = np.column_stack([self.x, self.y])

        return dict_

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################
    def to_json(self, path: str | pathlib.Path, encoding: str = "utf-8", **kwargs):
        import json

        kwargs = kwargs or dict()
        with open(path, 'w', encoding=encoding) as file:
            json.dump(self.to_dict(data_as_list=True), file, **kwargs)

    def to_csv(self, path: str | pathlib.Path, headers: bool = False, encoding: str = "utf-8"):
        kwargs = {"encoding": encoding}
        if headers:
            kwargs["headers"] = [self.x_label, self.y_label]
        np.savetxt(path, np.column_stack((self.x, self.y)), delimiter=",", **kwargs)

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        np.save(path, np.column_stack((self.x, self.y)), **kwargs)

    def to_npz(self, path: str | pathlib.Path, **kwargs):
        np.savez(path, x=self.x, y=self.y, name=self.name, x_label=self.x_label, y_label=self.y_label, **kwargs)

    def to_feather(self, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import numpy_to_feather

        headers = [self.x_label, self.y_label]
        numpy_to_feather(np.column_stack((self.x, self.y)), path, headers=headers)

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)
        if path.suffix == ".npz":
            return cls.from_npz(path)
        elif path.suffix == ".feather":
            return cls.from_feather(path)
        elif path.suffix == ".csv":
            return cls.from_csv(path)
        elif path.suffix == ".npy":
            return cls.from_npy(path)
        else:
            raise ValueError(f"Unsupported file type: {path.suffix}")

    @classmethod
    def from_json(cls, path: str | pathlib.Path, encoding: str = "utf-8", **kwargs):
        if isinstance(path, str):
            path = pathlib.Path(path)

        import json
        with open(path, 'r', encoding=encoding) as file:
            data = json.load(file)

        signal = np.array(data['y'])
        x = signal[:, 0]
        y = signal[:, 1]
        return cls(x, y, x_label=data['x_label'], y_label=data['y_label'], id_=data['id_'])

    @classmethod
    def from_csv(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        x, y, x_label, y_label = _load_csv(path)
        return cls(x, y, x_label=x_label, y_label=y_label)

    @classmethod
    def from_npy(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        x, y = np.load(str(path))
        return cls(x, y)

    @classmethod
    def from_npz(cls, path: str | pathlib.Path, **kwargs):
        npzfile = np.load(str(path), **kwargs)
        x, y, name, x_label, y_label = npzfile['x'], npzfile['y'], npzfile['name'], npzfile['x_label'], npzfile['y_label']
        return cls(x, y, x_label=x_label, y_label=y_label, name=name)

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


def _load_csv(path: pathlib) -> tuple[np.ndarray, np.ndarray, str | None, str | None]:
    import csv

    data = []
    x_label = None
    y_label = None

    with open(path, 'r') as file:
        csv_reader = csv.reader(file)

        first_row = next(csv_reader, None)
        if len(first_row) != 2:
            raise ValueError("Data not correct format. It should be in [n by 2] format.")
        if any(not is_number(cell) for cell in first_row):
            # If the first row contains non-numeric values, consider it as column labels
            x_label, y_label = first_row
        else:
            # If the first row contains numbers, treat them as data and set labels to None
            data.append([to_number(cell) for cell in first_row])

        for row in csv_reader:
            data.append([to_number(cell) for cell in row])

    data_array = np.array(data)
    return data_array[:, 0], data_array[:, 1], x_label, y_label


def is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def to_number(s: str) -> float | int:
    num = float(s)
    if int(num) == num:
        return int(num)
    else:
        return num
