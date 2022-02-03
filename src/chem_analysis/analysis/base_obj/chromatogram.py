from typing import Union

import numpy as np
import pandas as pd

from src.chem_analysis.analysis.base_obj.signal import Signal


class Chromatogram:
    """
    A grouping of Signals.
    """
    def __init__(self, data: Union[pd.DataFrame, Signal, list[Signal]]):
        super().__init__()
        if isinstance(data, pd.DataFrame):
            data = self._load_from_df(data)
        elif isinstance(data, Signal):
            data = {data.y_label if data.y_label is not None else "y_axis": data}
        elif isinstance(data, list):
            if all(isinstance(dat, Signal) for dat in data):
                data = {dat.y_label if dat.y_label is not None else f"y_axis{i}": dat for i, dat in enumerate(data)}
            else:
                raise ValueError("Invalid type in list")

        self.signals = data
        for k, v in self.signals.items():
            setattr(self, k, v)

    @property
    def labels(self):
        return list(self.signals.keys())

    @property
    def num_traces(self):
        return len(self.signals)

    def _load_from_df(self, df: pd.DataFrame) -> list[Signal]:
        """ Converts pandas dataframe into traces. """
        signals = []
        for col in df.columns:
            signals[col] = Signal(raw=df[col])

        return signals

    def to_dataframe(self):
        pass

    def baseline(self):
        pass

    def despike(self):
        pass

    def smooth(self):
        pass

    def auto_peak_picking(self):
        pass

    def stats(self):
        pass

    def plot(self):
        pass

    def plot_sep_y(self):
        pass


def local_run():
    from scipy.stats import norm
    nx = 1000
    ny = 3
    rv = norm(loc=nx/2, scale=10)
    x = np.linspace(0, nx, nx)
    y = np.empty((nx, ny))
    for i in range(ny):
        y[i] = np.linspace(0, nx, nx) + 20 * np.random.random(nx) + 5000 * rv.pdf(x) * np.random.random(1)

    df = pd.DataFrame(data=y, index=x)
    chro = Chromatogram(df)
    chro.plot()
    print("done")


if __name__ == "__main__":
    local_run()

    # path = r"C:\Users\nicep\Desktop\Reseach_Post\Data\Polyester\DW1_3\SEC\DW1-3-2[DW1-3].csv"
    # df = pd.read_csv(path, header=0, index_col=0)
    # df = df.iloc[:, :10]
    # df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    # df.index.names = ["time (min)"]