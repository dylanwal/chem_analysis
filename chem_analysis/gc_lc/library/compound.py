from __future__ import annotations
from typing import Iterable
import base64
import logging
from collections import OrderedDict
import ast

import numpy as np
import bigsmiles

from chem_analysis.config import global_config
from chem_analysis.mass_spec.ms_signal import MSSignal

logger = logging.getLogger(__name__)


OPTIMIZATION_KEY = {
    "name": "name",
    "datetime_created": "ctime",
    "datetime_updated": "utime",
    "compounds": "comps",
    "label": "label",
    "responses": "resps",
    "method": "method",
    "retention_time": "rtime",
    "response": "resp",
    "mass_spectrum": "ms",
    "notes": "notes",
    "mass_spectrum_encoding": "ms_e",
    "smiles": "smiles",
    "groups": "group",
    "cas": "cas",
    "density": "den",
    "boiling_temperature": "btemp",
    "parent": "parent"

}
OPTIMIZATION_KEY_REVERSE = {v: k for k, v in OPTIMIZATION_KEY.items()}


def is_optimized(keys: Iterable) -> bool:
    normal_count = 0
    optimized_count = 0
    for k in keys:
        if k in OPTIMIZATION_KEY:
            normal_count += 1
        if k in OPTIMIZATION_KEY_REVERSE:
            optimized_count += 1
    return optimized_count > normal_count


class NumpyEncoding:
    TEXT_ENCODING = "ASCII"
    DEFAULT_BYTES = "Base64"

    def __init__(self,
                 bytes_: str,
                 dtype: str,
                 shape: tuple[int, ...]
                 ):
        self.bytes_ = bytes_
        self.dtype = dtype
        self.shape = shape

    def __str__(self):
        return f"{self.bytes_};{self.dtype};{self.shape}"

    @classmethod
    def from_array(cls, array: np.ndarray) -> NumpyEncoding:
        return cls(cls.DEFAULT_BYTES, array.dtype.str, array.shape)

    @classmethod
    def from_string(cls, text: str) -> NumpyEncoding:
        text = text.split(";")
        return cls(bytes_=text[0], dtype=text[1], shape=ast.literal_eval(text[2]))


class CompoundResponse:
    # Update version in GCLibrary if changes made
    __slots__ = "method", "retention_time", "response", "mass_spectrum", "notes"
    RETENTION_TIME_UNIT = "min"

    def __init__(self,
                 method: str,
                 retention_time: int | float = None,  # min
                 response: int | float | None = None,
                 mass_spectrum: np.ndarray = None,
                 notes: str = None
                 ):
        self.method = method
        self.retention_time = retention_time
        self.response = response
        self.mass_spectrum = mass_spectrum
        self.notes = notes

    def __str__(self):
        text = f"{self.method}"
        text += f" ({self.retention_time} {self.RETENTION_TIME_UNIT})" or ""
        text += f" | ms_peaks: {self.mass_spectrum.shape[0]}" if self.mass_spectrum is not None else ""
        return text

    def __repr__(self):
        return self.__str__()

    def to_dict(self, sanitize: bool = False, binary: bool = False, optimize: bool = False) -> OrderedDict:
        dict_ = OrderedDict()
        vars_ = self.__slots__
        for k in vars_:
            attr = getattr(self, k)

            if sanitize:
                if k == 'mass_spectrum' and attr is not None:
                    if binary:
                        encoding = NumpyEncoding.from_array(attr)
                        dict_["mass_spectrum_encoding"] = str(encoding)
                        attr = base64.b64encode(attr.tobytes()).decode(NumpyEncoding.TEXT_ENCODING)
                    else:
                        attr = attr.tolist()

            if optimize:
                if k == "retention_time" and attr is not None and isinstance(attr, float):
                    attr = round(attr, global_config.sig_fig)
                if k == "response" and attr is not None and isinstance(attr, float):
                    attr = round(attr, global_config.sig_fig)
                if "mass_spectrum_encoding" in dict_:
                    dict_[OPTIMIZATION_KEY["mass_spectrum_encoding"]] = dict_.pop("mass_spectrum_encoding")

                k = OPTIMIZATION_KEY[k]

            dict_[k] = attr

        return dict_

    @classmethod
    def from_dict(cls, dict_: dict, optimized: bool = None) -> CompoundResponse:
        if optimized is None:
            optimized = is_optimized(dict_.keys())
        if optimized:
            dict_ = {OPTIMIZATION_KEY_REVERSE.get(k, k): v for k, v in dict_.items()}

        if dict_["mass_spectrum"] is not None:
            if isinstance(dict_["mass_spectrum"], str):
                numpy_encoding = NumpyEncoding.from_string(dict_.pop("mass_spectrum_encoding"))
                dict_["mass_spectrum"] = np.frombuffer(
                    base64.b64decode(dict_["mass_spectrum"]),
                    dtype=numpy_encoding.dtype,
                ).reshape(numpy_encoding.shape)
            elif isinstance(dict_["mass_spectrum"], np.ndarray):
                dict_["mass_spectrum"] = np.array(dict_["mass_spectrum"])

        return cls(**dict_)


class Compound:
    # Update version in GCLibrary if changes made
    __slots__ = "label", "responses", "groups", "smiles", "name", "cas", "density", "boiling_temperature", "parent"

    def __init__(self,
                 label: str,
                 responses: list[CompoundResponse],
                 groups: str | Iterable[str] = None,
                 smiles: str | None = None,
                 name: str | None = None,
                 cas: str | None = None,
                 density: int | float | None = None,
                 boiling_temperature: int | float | None = None,
                 parent: str = None
                 ):
        """

        Parameters
        ----------
        label:
            main id
        responses:

        groups:
            keywords that can be used to group compounds
        smiles:
            smile chemical representation
        name:
            name of chemical
        cas:
            cas number
        parent:
            if compound is from derivatization
        """
        self.label = label
        self.responses = responses
        if groups is not None and isinstance(groups, str):
            groups = [groups]
        self.groups = groups or []
        self.name = name or label
        self.cas = cas
        if isinstance(smiles, str):
            smiles = bigsmiles.BigSMILES(smiles)
        self.smiles = smiles
        self.density = density
        self.boiling_temperature = boiling_temperature
        self.parent = parent

    def __str__(self):
        text = f"{self.label}"
        text += f", {self.name}" or ""
        text += f" {self.groups}" or ""
        text += f" | responses: {len(self.responses)}"
        return text

    def __repr__(self):
        text = self.__str__()
        text += f", {self.cas}" or ""
        text += f", {self.smiles}" or ""
        return text

    @property
    def methods(self) -> list[str]:
        return list(response.method for response in self.responses)

    def get_ms(self) -> None | MSSignal:
        for response in self.responses:
            if response.mass_spectrum is not None:
                return MSSignal(x_raw=response.mass_spectrum[:, 0], data_raw=response.mass_spectrum[:, 1])

        return None

    def to_dict(self, sanitize: bool = False, binary: bool = False, optimize: bool = False) -> OrderedDict:
        dict_ = OrderedDict()
        vars_ = self.__slots__
        for k in vars_:
            attr = getattr(self, k)

            if sanitize:
                if k == 'smiles' and attr is not None:
                    attr = str(attr)
                if k == 'responses':
                    attr = [response.to_dict(sanitize=sanitize, binary=binary, optimize=optimize) for response in attr]

            if optimize:
                if k == "density" and attr is not None and isinstance(attr, float):
                    attr = round(attr, global_config.sig_fig)
                if k == "boiling_temperature" and attr is not None and isinstance(attr, float):
                    attr = round(attr, global_config.sig_fig)

                k = OPTIMIZATION_KEY[k]

            dict_[k] = attr

        return dict_

    @classmethod
    def from_dict(cls, dict_: dict, optimized: bool = None) -> Compound:
        if optimized is None:
            optimized = is_optimized(dict_.keys())
        if optimized:
            dict_ = {OPTIMIZATION_KEY_REVERSE.get(k, k): v for k, v in dict_.items()}

        if dict_["smiles"] is not None:
            dict_["smiles"] = bigsmiles.BigSMILES(dict_["smiles"])
        if dict_["responses"] is not None:
            dict_["responses"] = [CompoundResponse.from_dict(response) for response in dict_["responses"]]

        return cls(**dict_)
