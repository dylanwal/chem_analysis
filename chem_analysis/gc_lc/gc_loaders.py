import pathlib

from chem_analysis.gc_lc.gc_parameters import GCParameters
from chem_analysis.gc_lc.gc_signal import GCSignal
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
            ms_data_2d = MSSignal2D.from_list(data=ms_data.pop('data'), y=ms_data.pop('time'))
            gc_signal = GCMSSignal(ms_raw=ms_data_2d)
            gc_signal.parameters = GCParameters(**(ini_data | ms_data))

        if fid_data is not None:
            fid_signal = GCSignal(x_raw=fid_data.pop('time'), data_raw=fid_data.pop('data'))
            fid_signal.parameters = GCParameters(**(ini_data | fid_data))

        return gc_signal, fid_signal
