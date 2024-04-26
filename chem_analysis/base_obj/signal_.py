from typing import Sequence
import pathlib

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.processing.base import Processor
from chem_analysis.analysis.peak import PeakBounded


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
    _peak_type = PeakBounded

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

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        if path.suffix == ".csv":
            x, y, x_label, y_label = load_csv(path)
        elif path.suffix == ".feather":
            from chem_analysis.utils.feather_format import feather_to_numpy
            data, headers = feather_to_numpy(path)
            x, y = data[:, 0], data[:, 1]
            if headers[0] != "0":
                x_label = headers[0]
                y_label = headers[1]
            else:
                x_label = y_label = None
        elif path.suffix == ".npy":
            x, y = np.load(str(path))
            x_label = y_label = None
        else:
            raise NotImplementedError("File type currently not supported.")

        return cls(x, y, x_label=x_label, y_label=y_label)

    def to_feather(self, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import numpy_to_feather

        headers = [self.x_label, self.y_label]
        numpy_to_feather(np.column_stack((self.x, self.y)), path, headers=headers)

    def to_csv(self, path: str | pathlib.Path, headers: bool = False, encoding: str = "utf-8"):
        kwargs = {"encoding": encoding}
        if headers:
            kwargs["headers"] = [self.x_label, self.y_label]
        np.savetxt(path, np.column_stack((self.x, self.y)), delimiter=",", **kwargs)

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        np.save(path, np.column_stack((self.x, self.y)), **kwargs)


def load_csv(path: pathlib) -> tuple[np.ndarray, np.ndarray, str | None, str | None]:
    import csv

    data = []
    x_label = None
    y_label = None

    with open(path, 'r') as file:
        csv_reader = csv.reader(file)

        # Check if the first row contains numbers
        first_row = next(csv_reader, None)

        if len(first_row) != 2:
            raise ValueError("Data not correct format.")

        if any(cell.isalpha() for cell in first_row):
            # If the first row contains non-numeric values, consider it as column labels
            x_label, y_label = first_row
        else:
            # If the first row contains numbers, treat them as data and set labels to None
            data.append([float(cell) for cell in first_row])

        # Read the remaining rows
        for row in csv_reader:
            data.append([float(cell) for cell in row])

    # Convert data to NumPy array
    data_array = np.array(data)

    return data_array[:, 0], data_array[:, 1], x_label, y_label
