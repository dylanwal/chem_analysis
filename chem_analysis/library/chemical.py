from typing import Any
from collections import OrderedDict

from chem_analysis.library.attribute import Attribute
from chem_analysis.library.identifiers import Identifier


class Chemical:
    __slots__ = "name", "identifiers", "attributes", "id_"

    def __init__(self,
                 name: str,
                 identifiers: list[Identifier] | None = None,
                 attributes: list[Attribute] | None = None,
                 id_: int | None = None,
                 ):
        """

        Parameters
        ----------
        name:
        identifiers:
        attributes:
        """
        self.name = name
        self.attributes = attributes or []
        self.identifiers = identifiers or []
        self.id_ = id_

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

    def to_json(self, /, **kwargs) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["name"] = self.name
        dict_["id_"] = self.id_
        dict_["attributes"] = [attr.to_json(**kwargs) for attr in self.attributes]
        dict_["identifiers"] = [iden.to_json(**kwargs) for iden in self.identifiers]

        return dict_

    @classmethod
    def _from_JSON(cls, dict_: dict, /, **kwargs):
        dict_["attributes"] = [Attribute._from_JSON(attr, **kwargs) for attr in dict_["attributes"]]
        dict_["identifiers"] = [Identifier._from_JSON(iden, **kwargs) for iden in dict_["identifiers"]]
        return cls(**dict_)

    def match_identifier(self,
                         label: str,
                         type_: type[Identifier],
                         ) -> bool:
        idens = [iden for iden in self.identifiers if isinstance(iden, type_)]
        if len(idens) == 0:
            return False

        for iden in idens:
            if label in iden:
                return True

        return False

    def get_identifier(self, type_: type[Identifier], class_: str | None = None) -> Identifier:
        for iden in self.identifiers:
            if isinstance(iden, type_):
                if hasattr(iden, class_):
                    if class_ is not None and class_ == iden.class_:
                        return iden
                else:
                    return iden

    def get_attribute(self, type_: type[Attribute], class_: str | None = None) -> Attribute:
        for attr in self.attributes:
            if isinstance(attr, type_):
                if hasattr(attr, class_):
                    if class_ is not None and class_ == attr.class_:
                        return attr
                else:
                    return attr
