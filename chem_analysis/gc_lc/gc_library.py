from __future__ import annotations
import base64
import pathlib

import numpy as np
import bigsmiles

MAX_mz = 1000


class Compound:
    def __init__(self,
                 label: str,
                 group: str = None,
                 retention_time: int | float = None,
                 smiles: str | None = None,
                 name: str | None = None,
                 cas: str | None = None,
                 response: int | float | None = None,
                 mass_spectrum: np.ndarray = None,
                 ):
        self.label = label
        self.group = group
        self.name = name or label
        self.cas = cas
        if isinstance(smiles, str):
            smiles = bigsmiles.BigSMILES(smiles)
        self.smiles = smiles
        self.retention_time = retention_time
        self.response = response
        self.mass_spectrum = mass_spectrum

    def __str__(self):
        return f"{self.name}, {self.group}, {self.retention_time} min, {self.response}"

    def to_dict(self, sanitize: bool = False) -> dict:
        dict_ = {k: getattr(self, k) for k in vars(self) if not k.startswith("_")}

        if sanitize:
            if dict_['smiles'] is not None:
                dict_['smiles'] = str(dict_['smiles'])
            if dict_["mass_spectrum"] is not None:
                filtered_data = dict_["mass_spectrum"][dict_["mass_spectrum"] != 0]
                mz = np.nonzero(dict_["mass_spectrum"])
                dict_["mass_spectrum"] = np.column_stack([mz, filtered_data]).tolist()
                # dict_["mass_spectrum"] = base64.b64encode(obj.tobytes()).decode('utf-8')

        return dict_

    @classmethod
    def from_dict(cls, dict_: dict) -> Compound:
        if dict_["mass_spectrum"] is not None:
            data = np.array(dict_["mass_spectrum"])
            array_ = np.zeros(MAX_mz)
            index = (data[:, 0]).as_type('uint64')
            array_[index] = data[:, 1]
            dict_["mass_spectrum"] = array_
            # dict_["mass_spectrum"] = np.frombuffer(base64.b64decode(dict_["mass_spectrum"]), dtype=np.float64)
        if dict_["smiles"] is not None:
            dict_["smiles"] = bigsmiles.BigSMILES(dict_["smiles"])

        return cls(**dict_)


#######################################################################################################################
class GCLibrary:
    def __init__(self, compounds: list[Compound] = None, name: str = None, gc_method: str = None):
        self.name = name
        self.gc_method = gc_method
        self.compounds = []
        if compounds:
            for compound in compounds:
                self.add_compound(compound)

        # cache
        self._groups = None
        self._index = 0

    def __str__(self):
        return f"{len(self.compounds)} compounds"

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.compounds):
            compound = self.compounds[self._index]
            self._index += 1
            return compound
        else:
            raise StopIteration

    def _reset_cache(self):
        self._groups = None
        self._index = 0

    @property
    def groups(self) -> dict[str, list[Compound]]:
        if self._groups is None:
            groups = {}
            for compound in self.compounds:
                if compound.group not in groups:
                    groups[compound.group] = [compound]
                else:
                    groups[compound.group].append(compound)
            self._groups = groups

        return self._groups

    def add_compound(self, compound: Compound):
        self.compounds.append(compound)
        self._reset_cache()

    def find_by_name(self, name: str) -> Compound | None:
        for compound in self.compounds:  # TODO: add fuzzy matching
            if compound.name == name:
                return compound

    def find_by_label(self, label: str) -> Compound | None:
        for compound in self.compounds:
            if compound.label == label:
                return compound

    def find_by_cas(self, cas: str) -> Compound | None:
        for compound in self.compounds:
            if compound.cas == cas:
                return compound

    def find_by_smiles(self, smiles: str) -> Compound | None:
        # just text match # TODO: make SMART Search
        for compound in self.compounds:
            if compound.smiles == smiles:
                return compound

    # def get_n_nearest_compounds(self, retention_time: float, n: int = 1) -> list[Compound]:
    #     distance = []
    #     for i, peak in enumerate(self.peaks):
    #         distance[i] = abs(peak.retention_time - retention_time)
    #
    #     sort_index = np.argsort(distance)
    #     compounds = []
    #     for i in range(n):
    #         compounds.append(self.compounds[sort_index[i]])
    #     return compounds

    def to_dict(self, sanitize: bool = False) -> dict:
        dict_ = {k: getattr(self, k) for k in vars(self) if not k.startswith("_")}

        if sanitize:
            dict_['compounds'] = [comp.to_dict(sanitize) for comp in dict_['compounds']]

        return dict_

    def to_JSON(self, file_path: str | pathlib.Path):
        import json
        lib_dict = self.to_dict(sanitize=True)

        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)
        if file_path.suffix != ".json":
            file_path = file_path.with_suffix(".json")

        with open(file_path, 'w', encoding='UTF-8') as file:
            json.dump(lib_dict, file, indent=4)

    @classmethod
    def from_dict(cls, dict_: dict) -> GCLibrary:
        if dict_["compounds"] is not None:
            dict_["compounds"] = [Compound.from_dict(compound) for compound in dict_["compounds"]]
        return cls(**dict_)

    @classmethod
    def from_JSON(cls, file_path: str) -> GCLibrary:
        import json

        with open(file_path, 'r', encoding='UTF-8') as file:
            lib = cls.from_dict(json.load(file))
        return lib

    def to_picking_library(self):
        from chem_analysis.analysis.peak_picking.library_search import PickingLibrary, PeakForPicking
        peaks = []
        for compound in self:
            if compound.retention_time is not None:
                peaks.append(PeakForPicking(compound.retention_time))

        return PickingLibrary(peaks)
