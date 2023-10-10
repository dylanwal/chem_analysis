class NMR:
    def __init__(self, parameters: NMRParameters = None):
        self._fid = None
        self._spectrum = None
        self._FID_time_axis = None
        self._ppm_axis = None
        self.parameters = parameters
        self.processing = Processor()

    def __repr__(self):
        return f"{self.parameters.nucleus} (in {self.parameters.solvent})"

    def __str__(self):
        return self.__repr__()

    @property
    def FID(self) -> np.ndarray:
        return self._fid

    @property
    def FID_real(self) -> np.ndarray:
        """ y axis of FID for visualization """
        return np.real(self._fid)

    @property
    def FID_time_axis(self) -> np.ndarray:
        """ x axis of FID """
        if self._FID_time_axis is None:
            self._FID_time_axis = np.linspace(
                0,
                (self.parameters.sizeTD2 - 1) * self.parameters.dwell_time,
                self.parameters.sizeTD2
            )

        return self._FID_time_axis

    @property
    def spectrum(self) -> np.ndarray:
        return 1

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
    def from_bruker(cls, path: pathlib.Path) -> NMR:
        # load data from file
        parameters = parse_acqus_file(path / "acqus")
        data = get_fid(path / "fid")

        # construct NMR object
        nmr = NMR(parameters=parameters)
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


def get_fid(path: pathlib.Path, endianess: str = "<", dtype: np.dtype = np.dtype("f8")) -> np.ndarray:
    dtype_ = dtype.newbyteorder(endianess)
    with open(path, mode='rb') as f:
        d = f.read()
        return np.frombuffer(d, dtype=dtype_)
