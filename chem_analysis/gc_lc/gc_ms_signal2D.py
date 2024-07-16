import logging
import pathlib
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.mass_spec.ms_signal_2D import MSSignal3D
from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal

logger = logging.getLogger(__name__)


class GCMSSignal2D(Signal2D):
    _signal = GCMSSignal

    def __init__(self,
                 ms_raw: MSSignal3D,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "time"
        z_label = z_label or "intensity"
        z_raw = np.sum(ms_raw, axis=2)
        super().__init__(ms_raw.x_raw, ms_raw.y_raw, z_raw, x_label, y_label, z_label, name)
        self.ms_raw = ms_raw

    def get_signal(self, y_index: int, processed: bool = False) -> GCMSSignal:
        return super().get_signal(y_index, processed)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[GCMSSignal],
                     y: np.ndarray = None,
                     y_label: str = None,
                     min_ms: int = None,
                     max_ms: int = None,
                     ):  # -> Signal2D
        if y is not None and len(y.shape) != 1 and y.shape[0] == len(signals):
            raise ValueError("The number of signals must be the same as the number of y points.\n"
                             f"\tnumber of signals: {len(signals)}\n\tnumber of y points:{y.shape[0]}")
        if y is None:
            y = np.arange(len(signals))

        x_label = signals[0].x_label
        z_label = signals[0].y_label
        # for i, sig in enumerate(signals):
        #     if len(sig.x) != len(x):
        #         raise ValueError(f"Signal {i} has different number of data points that Signal 1.\n")
        #     if not np.all(np.isclose(sig.x, x, rtol=0.01)):
        #         raise ValueError(f"Signal {i} has a different x-axis than first signal.")

        ms_raw = MSSignal3D.from_signals([sig.ms_raw for sig in signals], y, min_ms, max_ms)
        return cls(ms_raw=ms_raw, x_label=x_label, y_label=y_label, z_label=z_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################
    #
    # @classmethod
    # def from_file(cls, path: str | pathlib.Path):
    #     if isinstance(path, str):
    #         path = pathlib.Path(path)
    #
    #     if path.suffix == ".csv":
    #         raise NotImplementedError()
    #     elif path.suffix == ".feather":
    #         raise NotImplementedError()
    #     elif path.suffix == ".npz":
    #         npzfile = np.load(str(path))
    #         x, y, ms_raw = npzfile['x'], npzfile['y'], npzfile['ms_raw']
    #         x_label = y_label = z_label = None
    #     else:
    #         raise NotImplemented("File type currently not supported.")
    #
    #     return cls(x_raw=x, y_raw=y, ms_raw=ms_raw, x_label=x_label, y_label=y_label, z_label=z_label)
    #
    # def to_feather(self, path: str | pathlib.Path):
    #     raise NotImplementedError()
    #
    # def to_csv(self, path: str | pathlib.Path, **kwargs):
    #     raise NotImplementedError()
    #
    # def to_npy(self, path: str | pathlib.Path, **kwargs):
    #     import warnings
    #     warnings.warn("Redirected 'GCMSSignal2D.to_npy()' to 'GCMSSignal2D.to_npz()'.")
    #     self.to_npz(path)
    #
    # def to_npz(self, path: str | pathlib.Path, **kwargs):
    #     """Save an array to a binary file in NumPy ``.npz`` format."""
    #     if isinstance(path, str):
    #         path = pathlib.Path(path)
    #     if path.suffix != ".npz":
    #         path = path.with_suffix(".npz")
    #     np.savez(path, ms_raw=self.ms_raw, **kwargs)
