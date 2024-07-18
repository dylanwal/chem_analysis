import pathlib

import numpy as np

from chem_analysis.gc_lc.gc_parameters import GCParameters
from chem_analysis.gc_lc.gc_signal import GCSignal
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D
from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal


class GCParser:
    @classmethod
    def from_Agilent_D_folder(cls, folder_path: str | pathlib.Path) -> tuple[GCMSSignal | None, GCSignal | None]:
        if not isinstance(folder_path, pathlib.Path):
            folder_path = pathlib.Path(folder_path)

        from chem_analysis.gc_lc.parsers.agilent_folder import parse_D_folder

        ini_data, ms_data, fid_data = parse_D_folder(folder_path)

        gc_signal, fid_signal = None, None
        if ms_data is not None:
            data = ms_data.pop('data')
            for i, sig in enumerate(data):
                if len(sig[0]) == 0:
                    data[i] = [np.array([0]), np.array([0])]

            ms_list = [MSSignal(x=ms[0], y=ms[1]) for ms in data]
            ms_data_2d = MSSignal2D.from_signals(ms_list, y=ms_data.pop('time'))
            ms_data_2d.name = ms_data.pop('sample_name')
            gc_signal = GCMSSignal(ms=ms_data_2d)
            gc_signal.parameters = GCParameters(**(ini_data | ms_data))

        if fid_data is not None:
            fid_signal = GCSignal(
                x=fid_data.pop('time'),
                y=fid_data.pop('data'),
                name=fid_data.pop('sample_name')
            )
            fid_signal.parameters = GCParameters(**(ini_data | fid_data))

        return gc_signal, fid_signal
