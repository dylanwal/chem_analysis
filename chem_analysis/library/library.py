from __future__ import annotations

import datetime
import logging
import pathlib
from collections import OrderedDict
from typing import Any

from chem_analysis.library.chemical import Chemical
from chem_analysis.library.identifiers import Identifier

logger = logging.getLogger(__name__)


class Library:
    VERSION = 1  # update if anything changes here or in Chemical

    def __init__(self,
                 name: str,
                 chemicals: list[Chemical] = None,
                 datetime_created: datetime.datetime | str | None = None,
                 datetime_modified: datetime.datetime | str | None = None
                 ):
        self.name = name
        if isinstance(datetime_created, str):
            datetime_created = datetime.datetime.fromisoformat(datetime_created)
        if isinstance(datetime_modified, str):
            datetime_created = datetime.datetime.fromisoformat(datetime_modified)
        self.datetime_created = datetime_created or datetime.datetime.now()
        self.datetime_modified = datetime_modified or datetime.datetime.now()

        self._chemicals = []
        if chemicals:
            for chem in chemicals:
                self.add_chemical(chem)

    def __str__(self):
        return f"{len(self.chemicals)} chemicals"

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.chemicals):
            chemicals = self.chemicals[self._index]
            self._index += 1
            return chemicals
        else:
            raise StopIteration

    def __contains__(self, item: Chemical):
        if item in self.chemicals:
            return True
        return False

    @property
    def chemicals(self) -> list[Chemical]:
        return self._chemicals

    def add_chemical(self, chemicals: Chemical):
        if chemicals in self:
            logging.warning(f"Can't add duplicate chemicals. The original chemicals in {self} is retained. "
                            f"\n Chemical: {chemicals}")

        self.chemicals.append(chemicals)

    def delete_chemical(self, chemical: Chemical | str | int):
        index = None
        if isinstance(chemical, str):
            chemical = self.find_by_name(chemical)
        if isinstance(chemical, Chemical):
            for i, chem in enumerate(self):
                if chemical is chem:
                    index = i
                    break

        if index is None:
            raise ValueError(f"'chemicals' not found in library.\n\tchemicals: {str(chemical)}")

        self._chemicals.pop(index)

    def add_library(self, lib: Library):
        for comp in lib:
            self.add_chemical(comp)

    def find_by_name(self, name: str) -> Chemical | None:
        for chemicals in self.chemicals:
            if chemicals.name == name:
                return chemicals

    def find_by_identifier(self,
                           label: str,
                           class_: type[Identifier],
                           count: int = 1
                           ) \
            -> list[Chemical]:
        matches = []
        for i, chem in enumerate(self.chemicals):
            chem.match_identifier(label, class_)

            if len(matches) == count:
                logger.info(f"Reached count limit {count}. Searched {i}/{len(self.chemicals)}")
                break

        return matches

    def to_dict(self) -> OrderedDict:
        dict_ = OrderedDict()
        dict_["name"] = self.name
        dict_["datetime_created"] = self.datetime_created.isoformat()
        dict_["datetime_modified"] = self.datetime_modified.isoformat()
        dict_["chemicals"] = [chem.to_dict() for chem in self.chemicals]
        return dict_

    def to_JSON(self,
                file_path: str | pathlib.Path,
                *,
                numpy_encoding: str = "list",
                json_kwargs: dict[str, Any] = None,
                overwrite: bool = False
                ):
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)
        if file_path.suffix != ".json":
            file_path = file_path.with_suffix(".json")
        if file_path.exists() and not overwrite:
            raise ValueError("Library file already exists. Set 'overwrite' to true.")

        import json
        lib_dict = OrderedDict()
        lib_dict["version"] = self.VERSION
        lib_dict["encoding"] = "UTF-8"
        lib_dict["encoding_numpy"] = numpy_encoding
        lib_dict["name"] = self.name
        lib_dict["datetime_created"] = self.datetime_created.isoformat()
        lib_dict["datetime_updated"] = datetime.datetime.now().isoformat()

        chems = []
        for chem in self.chemicals:
            chems.append(chem.to_json(numpy_encoding))

        if json_kwargs is None:
            json_kwargs = {}
        if 'indent' not in json_kwargs:
            json_kwargs['indent'] = 2

        with open(file_path, 'w', encoding='UTF-8') as file:
            json.dump(lib_dict, file, **json_kwargs)

    def to_pickle(self, file_path: str | pathlib.Path, *, overwrite: bool = False):
        import pickle
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)
        if file_path.suffix != ".pkl":
            file_path = file_path.with_suffix(".pkl")
        if file_path.exists() and not overwrite:
            raise ValueError("Library file already exists. Set 'overwrite' to true.")

        with open(file_path, 'wb') as file:
            pickle.dump(self, file)

    @classmethod
    def from_pickle(cls, file_path: str | pathlib.Path) -> Library:
        import pickle
        with open(file_path, 'rb') as file:
            loaded_obj = pickle.load(file)
        return loaded_obj

    @classmethod
    def from_JSON(cls, file_path: str | pathlib.Path) -> Library:
        import json

        with open(file_path, 'r', encoding='UTF-8') as file:
            lib = json.load(file)

        version = lib.pop("version")
        if version != cls.VERSION:
            raise ValueError(f"Library version {version} does not match library version {cls.VERSION}")
        lib.pop("encoding")
        numpy_encoding = lib.pop("encoding_numpy")

        lib["compounds"] = [Chemical._from_JSON(chem, numpy_encoding) for chem in lib["chemical"]]

        return cls(**lib)
