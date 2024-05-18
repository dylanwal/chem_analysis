from __future__ import annotations

import datetime
import logging
import pathlib
from collections import OrderedDict, defaultdict
from typing import Any

from chem_analysis.gc_lc.library.compound import Compound

logger = logging.getLogger(__name__)


class GCLibrary:
    VERSION = 1  # update if anything changes here or in Compound

    def __init__(self,
                 name: str,
                 compounds: list[Compound] = None,
                 datetime_created: datetime.datetime = None,
                 datetime_updated: datetime.datetime = None
                 ):
        self.name = name
        self.datetime_created = datetime_created or datetime.datetime.now()
        self.datetime_updated = datetime_updated or datetime.datetime.now()  # updated in .to_JSON
        self._compounds = []
        if compounds:
            for compound in compounds:
                self.add_compound(compound)

        # cache values
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

    def __contains__(self, item: Compound):
        if item in self.compounds:
            return True
        return False

    @property
    def compounds(self) -> list[Compound]:
        return self._compounds

    def _reset_cache(self):
        self._groups = None
        self._index = 0

    def groups(self) -> dict[str, list[Compound]]:
        if self._groups is None:
            groups = defaultdict(list)
            for compound in self.compounds:
                for group in compound.groups:
                    groups[group].append(compound)
            self._groups = groups

        return self._groups

    def add_compound(self, compound: Compound):
        if compound in self:
            logging.warning(f"Can't add duplicate compound. The original compound in {self} is retained. "
                            f"\n Compound: {compound}")

        self.compounds.append(compound)
        self._reset_cache()

    def add_library(self, lib: GCLibrary):
        for comp in lib:
            self.add_compound(comp)

    def get_group(self, group: str) -> list[Compound]:
        return self._groups[group]

    def find_by_name(self, name: str) -> Compound | None:
        for compound in self.compounds:
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
            if str(compound.smiles) == smiles:
                return compound

    def to_dict(self, sanitize: bool = False, binary: bool = False) -> OrderedDict:
        dict_ = OrderedDict()
        vars_ = list(filter(lambda x: not x.startswith("_"), vars(self))) + ["_compounds"]
        for k in vars_:
            attr = getattr(self, k)
            if sanitize:
                if k == "_compounds":
                    k = "compounds"
                    attr = [comp.to_dict(sanitize, binary) for comp in attr]
                if isinstance(attr, datetime.datetime):
                    attr = attr.isoformat()

            dict_[k] = attr

        return dict_

    def to_JSON(self, file_path: str | pathlib.Path, binary: bool = False, json_kwargs: dict[str, Any] = None):
        import json
        lib_dict = self.to_dict(sanitize=True, binary=binary)
        lib_dict["datetime_updated"] = datetime.datetime.now().isoformat()
        lib_dict["version"] = self.VERSION

        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)
        if file_path.suffix != ".json":
            file_path = file_path.with_suffix(".json")

        if json_kwargs is None:
            json_kwargs = {}
        if 'indent' not in json_kwargs:
            json_kwargs['indent'] = 4

        with open(file_path, 'w', encoding='UTF-8') as file:
            json.dump(lib_dict, file, **json_kwargs)

    @classmethod
    def from_dict(cls, dict_: dict) -> GCLibrary:
        if dict_['version'] != cls.VERSION:
            raise ValueError(f"Version {dict_['version']} does not match the current version {cls.VERSION}")
        dict_.pop("version")

        if dict_["compounds"] is not None:
            dict_["compounds"] = [Compound.from_dict(compound) for compound in dict_["compounds"]]
        return cls(**dict_)

    @classmethod
    def from_JSON(cls, file_path: str | pathlib.Path) -> GCLibrary:
        import json

        with open(file_path, 'r', encoding='UTF-8') as file:
            lib = cls.from_dict(json.load(file))
        return lib
