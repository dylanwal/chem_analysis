import pathlib

import numpy as np

from chem_analysis.nmr.parameters import NMRParameters
from chem_analysis.processing.base import Processor
from chem_analysis.base_obj.signal_ import Signal


class NMRFID(Signal):
    """
    Free Induction Decay (FID)
    """
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = "time",
                 y_label: str = "signal",
                 name: str = None,
                 id_: int = None
                 ):
        super().__init__(x, y, x_label, y_label, name, id_)

    @property
    def FID_real(self) -> np.ndarray:
        """ y axis of FID for visualization """
        return np.real(self.y)


class NMRSignal(Signal):
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 parameters: NMRParameters = None,
                 x_label: str = "time",
                 y_label: str = "signal",
                 name: str = None,
                 id_: int = None
    ):
        super().__init__(x, y, x_label, y_label, name, id_)
        self._fid: NMRFID | None = None
        self.parameters = parameters
        self.processing = Processor()

    def __repr__(self):
        return f"{self.parameters.type_.name} (in {self.parameters.solvent})"

    def __str__(self):
        return self.__repr__()

    @property
    def ppm_axis(self) -> np.ndarray:
        if self._ppm_axis is None:
            self._ppm_axis = np.linspace(
                -self.parameters.spectral_width / 2,
                self.parameters.spectral_width / 2,
                self.parameters.sizeTD2
            )

        return self._ppm_axis

    def load_from_raw_FID_data(self, data: np.ndarray):
        _real = data[0:self.parameters.sizeTD2 * 2:2]
        _imag = np.multiply(data[1:self.parameters.sizeTD2 * 2 + 1:2], 1j)
        self._fid = np.add(_real, _imag)

    @classmethod
    def from_bruker(cls, path: pathlib.Path) -> NMRSignal:
        if isinstance(path, str):
            path = pathlib.Path(path)
        from chem_analysis.nmr.parse_bruker import parse_bruker_folder
        # load data from file
        fid_data, parameters = parse_brucker_folder

        # construct NMR object
        nmr = NMRSignal(parameters=parameters)
        nmr.load_from_raw_FID_data(data)
        return nmr

    @classmethod
    def from_spinsolve(cls, path: pathlib.Path | str) -> NMRSignal:
        if isinstance(path, str):
            path = pathlib.Path(path)
        from chem_analysis.nmr.parse_spinsolve import
        # load data from file
        parameters = parse_acqus_file(path / "acqus")
        data = get_fid(path / "fid")

        # construct NMR object
        nmr = NMRSignal(parameters=parameters)
        nmr.load_from_raw_FID_data(data)
        return nmr

    def default_processing(self):
        self.processing.operationStack = [
            OPS.LeftShift(self.parameters.shift_points),
            OPS.LineBroadening(0.0),
            OPS.FourierTransform(),
            OPS.Phase0D(-90),
            OPS.Phase1D(self.parameters.shift_points, unit="time")
        ]



