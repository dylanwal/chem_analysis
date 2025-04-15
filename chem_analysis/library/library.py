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

    def add_chemical(self, chemical: Chemical):
        if chemical in self:
            logging.warning(f"Can't add duplicate chemicals. The original chemicals in {self} is retained. "
                            f"\n Chemical: {chemical}")
        if chemical.id_ is None:
            chemical.id_ = len(self.chemicals)
        self.chemicals.append(chemical)

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
                           type_: type[Identifier],
                           class_: str | None = None,
                           count: int = 1
                           ) \
            -> list[Chemical]:
        matches = []
        for i, chem in enumerate(self.chemicals):
            if chem.match_identifier(label, type_, class_):
                matches.append(chem)

            if len(matches) == count:
                logger.info(f"Reached count limit {count} for search '{type_.__name__}.{label}'. Searched {i}/{len(self.chemicals)}")
                break

        return matches

    def to_dict(self) -> OrderedDict:
        dict_ = OrderedDict()
        dict_["name"] = self.name
        dict_["datetime_created"] = self.datetime_created.isoformat()
        dict_["datetime_modified"] = self.datetime_modified.isoformat()
        dict_["chemicals"] = [chem.to_dict() for chem in self.chemicals]
        return dict_

    def write_json(self,
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
        lib_dict["datetime_modified"] = datetime.datetime.now().isoformat()

        chems = []
        for chem in self.chemicals:
            chems.append(chem.write_json(numpy_encoding=numpy_encoding))

        lib_dict['chemicals'] = chems

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

    def to_dataframe(self):
        import polars as pl
        from collections import defaultdict, Counter

        def merge_df(dfs: list[pl.DataFrame]) -> pl.DataFrame:
            # Collect all column names and their observed types
            col_types = defaultdict(list)
            for df in dfs:
                for name, dtype in zip(df.columns, df.dtypes):
                    col_types[name].append(dtype)

            # Decide on the common dtype for each column (most frequent type wins)
            common_types = {}
            for col, types in col_types.items():
                most_common = Counter(types).most_common(1)[0][0]
                common_types[col] = most_common

            all_columns = set(common_types.keys())

            def pad_and_cast(df: pl.DataFrame) -> pl.DataFrame:
                # Add missing columns with correct dtype
                missing = [
                    pl.lit(None, dtype=common_types[col]).alias(col)
                    for col in all_columns if col not in df.columns
                ]
                df = df.with_columns(missing)

                # Cast existing columns to common types if needed
                for col in df.columns:
                    if df[col].dtype != common_types[col]:
                        df = df.with_columns([df[col].cast(common_types[col])])

                return df.select(sorted(all_columns))  # Optional: sort for consistent order

            dfs_padded = [pad_and_cast(df) for df in dfs]
            return pl.concat(dfs_padded, how="vertical")

        dfs = [chem.to_dataframe() for chem in self.chemicals]
        return merge_df(dfs)

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

        lib["chemicals"] = [Chemical._from_JSON(chem, numpy_encoding=numpy_encoding) for chem in lib["chemicals"]]

        return cls(**lib)
