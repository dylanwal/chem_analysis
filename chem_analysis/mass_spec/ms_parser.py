import pathlib


class MS_parser:
    @classmethod
    def from_Agilent_D_folder(cls, folder_path: str | pathlib.Path):
        if not isinstance(folder_path, pathlib.Path):
            folder_path = pathlib.Path(folder_path)

        from chem_analysis.mass_spec.parsers.agilent_folder import parse_D_folder

        ini_data, ms_data, fid_data = parse_D_folder(folder_path)
