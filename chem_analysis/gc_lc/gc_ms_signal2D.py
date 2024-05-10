import pathlib
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal


class GCMSSignal2D(Signal2D):
    _signal = GCMSSignal

    def __init__(self,
                 x_raw: np.ndarray,  # x
                 y_raw: np.ndarray,  # y
                 ms_raw: np.ndarray,  # [y, x, n]
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "time"
        z_label = z_label or "intensity"
        data_raw = np.sum(ms_raw, axis=2)
        super().__init__(x_raw, y_raw, data_raw, x_label, y_label, z_label, name)
        self.ms_raw = ms_raw

    def get_signal(self, y_index: int, processed: bool = False) -> GCMSSignal:
        if processed:
            sig = self._signal(x_raw=self.x, ms_raw=self.ms_raw[y_index], x_label=self.x_label,
                               y_label=self.y_label, name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
        else:
            sig = self._signal(x_raw=self.x_raw, ms_raw=self.ms_raw[y_index], x_label=self.x_label,
                               y_label=self.y_label, name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
            sig.processor = self.processor.get_copy()
        sig.y_value = self.y[y_index]
        return sig

    @classmethod
    def from_signals(cls, signals: Sequence[GCMSSignal], y: np.ndarray = None):  # -> Signal2D
        if y is not None and len(y.shape) != 1 and y.shape[0] == len(signals):
            raise ValueError("The number of signals must be the same as the number of y points.\n"
                             f"\tnumber of signals: {len(signals)}\n\tnumber of y points:{y.shape[0]}")
        if y is None:
            y = np.arang(len(signals))

        x = signals[0].x
        x_label = signals[0].x_label
        z_label = signals[0].y_label
        for i, sig in enumerate(signals):
            if not np.all(np.isclose(sig.x, x, rtol=0.01)):
                raise ValueError(f"Signal {i} has a different x-axis than first signal.")
        ms_raw = np.array([s.ms_raw for s in signals])
        return cls(x_raw=signals[0].x, y_raw=y, ms_raw=ms_raw, x_label=x_label, z_label=z_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)

        if path.suffix == ".csv":
            raise NotImplementedError()
        elif path.suffix == ".feather":
            raise NotImplementedError()
        elif path.suffix == ".npz":
            npzfile = np.load(str(path))
            x, y, ms_raw = npzfile['x'], npzfile['y'], npzfile['ms_raw']
            x_label = y_label = z_label = None
        else:
            raise NotImplemented("File type currently not supported.")

        return cls(x_raw=x, y_raw=y, ms_raw=ms_raw, x_label=x_label, y_label=y_label, z_label=z_label)

    def to_feather(self, path: str | pathlib.Path):
        raise NotImplementedError()

    def to_csv(self, path: str | pathlib.Path, **kwargs):
        raise NotImplementedError()

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        import warnings
        warnings.warn("Redirected 'GCMSSignal2D.to_npy()' to 'GCMSSignal2D.to_npz()'.")
        self.to_npz(path)

    def to_npz(self, path: str | pathlib.Path, **kwargs):
        """Save an array to a binary file in NumPy ``.npz`` format."""
        if isinstance(path, str):
            path = pathlib.Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")
        np.savez(path, x=self.x, y=self.y, ms_raw=self.ms_raw, **kwargs)
