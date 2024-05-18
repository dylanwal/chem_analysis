from __future__ import annotations
from typing import Iterable
import base64
import logging
from collections import OrderedDict

import numpy as np
import bigsmiles

logger = logging.getLogger(__name__)


class CompoundResponse:
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

    def to_dict(self, sanitize: bool = False, binary: bool = False) -> OrderedDict:
        dict_ = OrderedDict()
        vars_ = filter(lambda x: not x.startswith("_"), vars(self))
        for k in vars_:
            attr = getattr(self, k)

            if sanitize:
                if k == 'mass_spectrum' and attr is not None:
                    if binary:
                        attr = base64.b64encode(attr.tobytes()).decode('utf-8')
                    else:
                        attr = attr.tolist()

            dict_[k] = attr

        return dict_

    @classmethod
    def from_dict(cls, dict_: dict) -> CompoundResponse:
        if dict_["mass_spectrum"] is not None:
            if isinstance(dict_["mass_spectrum"], str):
                dict_["mass_spectrum"] = np.frombuffer(base64.b64decode(dict_["mass_spectrum"]))
            elif isinstance(dict_["mass_spectrum"], np.ndarray):
                dict_["mass_spectrum"] = np.array(dict_["mass_spectrum"])

        return cls(**dict_)


class Compound:
    __slots__ = "label", "responses", "groups", "smiles", "name", "cas"

    def __init__(self,
                 label: str,
                 responses: list[CompoundResponse],
                 groups: str | Iterable[str] = None,
                 smiles: str | None = None,
                 name: str | None = None,
                 cas: str | None = None,
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
        """
        self.label = label
        self.responses = responses
        if groups is not None and not isinstance(groups, str):
            groups = [groups]
        self.groups = groups or []
        self.name = name or label
        self.cas = cas
        if isinstance(smiles, str):
            smiles = bigsmiles.BigSMILES(smiles)
        self.smiles = smiles

    def __str__(self):
        text = f"{self.label}"
        text += f", {self.name}" or ""
        text += f"( {self.groups})" or ""
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

    def to_dict(self, sanitize: bool = False, binary: bool = False) -> OrderedDict:
        dict_ = OrderedDict()
        vars_ = filter(lambda x: not x.startswith("_"), vars(self))
        for k in vars_:
            attr = getattr(self, k)

            if sanitize:
                if k == 'smiles' and attr is not None:
                    attr = str(attr)
                if k == 'responses':
                    attr = [response.to_dict(sanitize, binary) for response in attr]

            dict_[k] = attr

        return dict_

    @classmethod
    def from_dict(cls, dict_: dict) -> Compound:
        if dict_["smiles"] is not None:
            dict_["smiles"] = bigsmiles.BigSMILES(dict_["smiles"])
        if dict_["responses"] is not None:
            dict_["responses"] = [CompoundResponse.from_dict(response) for response in dict_["responses"]]

        return cls(**dict_)
