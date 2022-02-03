from typing import Union, Tuple, List, Any
import re

import pandas as pd
import plotly.graph_objs as go

from src.chem_analysis.analysis.sec.trace import Trace





class Chromatogram:
    """
    traces need
    """
    def __init__(self, data: Union[pd.DataFrame, Trace, list[Trace]]):
        super().__init__()
        if isinstance(data, pd.DataFrame):
            data = self._load_from_df(data)
        elif isinstance(data, Trace):
            data = {data.y_label if data.y_label is not None else "y_axis": data}
        elif isinstance(data, list):
            if all(isinstance(dat, Trace) for dat in data):
                data = {dat.y_label if dat.y_label is not None else f"y_axis{i}": dat for i, dat in enumerate(data)}
            else:
                raise ValueError("Invalid type in list")

        self.traces = data
        for k, v in self.traces.items():
            setattr(self, k, v)

    @property
    def labels(self):
        return list(self.traces.keys())

    @property
    def num_traces(self):
        return len(self.traces)

    def _load_from_df(self, df: pd.DataFrame):
        """ Converts pandas dataframe into traces. """
        traces = {}
        x_label, x_unit = parse_label(df.index.name)
        for col in df.columns:
            y_label, y_unit = parse_label(col)
            traces[col] = Trace(
                x=df.index.to_numpy(),
                y=df[col].to_numpy(),
                x_label=x_label,
                y_label=y_label,
                x_unit=x_unit,
                y_unit=y_unit
            )

        return traces

    def to_dataframe(self, level=):





def local_run():
    path = r"C:\Users\nicep\Desktop\Reseach_Post\Data\Polyester\DW1_3\SEC\DW1-3-2[DW1-3].csv"
    df = pd.read_csv(path, header=0, index_col=0)
    df = df.iloc[:, :10]
    df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    df.index.names = ["time (min)"]

    chro = Chromatogram(df)
    chro.plot_sep_y()


if __name__ == "__main__":
    local_run()
