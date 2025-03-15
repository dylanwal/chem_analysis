from typing import Any
from collections import OrderedDict

from chem_analysis.library.attribute import Attribute
from chem_analysis.library.identifiers import Identifier


class Chemical:
    __slots__ = "name", "identifiers", "attributes", "id_", "datetime_created", "datetime_modified"

    def __init__(self,
                 name: str,
                 identifiers: list[Identifier] | None = None,
                 attributes: list[Attribute] | None = None,
                 ):
        """

        Parameters
        ----------
        name:
        identifiers:
        attributes:
        """
        self.name = name
        self.id_ = id
        self.attributes = attributes or None
        self.identifiers = identifiers or None

    def __str__(self):
        text = f"{self.name}"
        return text

    def __repr__(self):
        text = f"{self.name}"
        return text

    def to_dict(self) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["name"] = self.name
        dict_["id_"] = self.id_
        dict_["attributes"] = [attr.to_dict() for attr in self.attributes]
        dict_["identifiers"] = [iden.to_dict() for iden in self.identifiers]
        return dict_

    def to_json(self, numpy_encoding: str = 'list') -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["name"] = self.name
        dict_["id_"] = self.id_
        dict_["attributes"] = [attr.to_json(numpy_encoding) for attr in self.attributes]
        dict_["identifiers"] = [iden.to_json(numpy_encoding) for iden in self.identifiers]

        return dict_

    @classmethod
    def _from_JSON(cls, dict_: dict, numpy_encoding: str = 'list'):
        dict_["attributes"] = [Attribute._from_JSON(attr, numpy_encoding) for attr in dict_["attributes"]]
        dict_["identifiers"] = [Identifier._from_JSON(iden, numpy_encoding) for iden in dict_["identifiers"]]
        return cls(**dict_)

    def match_identifier(self,
                         label: str,
                         class_: type[Identifier],
                         ) -> bool:
        idens = [iden for iden in self.identifiers if isinstance(iden, class_)]
        if len(idens) == 0:
            return False

        for iden in idens:
            if label in iden:
                return True

        return False
