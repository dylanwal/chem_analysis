import abc
from typing import Any
from collections import OrderedDict
import inspect
import re

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = inspect.getmembers(class_, lambda a: not (inspect.isroutine(a)))
    list_of_attrs = [a for a in attrs if not (a[0].startswith('_'))]
    return {k: v for k, v in list_of_attrs}


class Identifier(MixinSubClassList, abc.ABC):
    __slots__ = ("value",)

    def __init__(self, value: Any):
        self.value = value

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return self.__str__()

    def __contains__(self, key: str) -> bool:
        if isinstance(self.value, str):
            return key == self.value
        if isinstance(self.value, list):
            for v in self.value:
                if key in v:
                    return True

        # try to work with custom types
        try:
            if key == self.value:
                return True
        except TypeError:
            pass

        return False

    def to_dict(self, remove_nones: bool = False) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["type"] = type(self).__name__
        dict_["value"] = self.value

        attrs = get_class_instance_attributes(self)
        attrs.pop('value')
        if remove_nones:
            attrs = {k: v for k, v in attrs.items() if v is not None}
        for k, v in attrs.items():
            dict_[k] = v
        return dict_

    def to_json(self, /, **kwargs) -> OrderedDict[str, Any]:
        return self.to_dict(remove_nones=True)

    @classmethod
    def _from_JSON(cls, dict_: OrderedDict[str, Any], /, **kwargs):
        class_ = dict_.pop("type")
        for k in cls.sub_classes():
            if k.__name__ == class_:
                class_ = k
        return class_._from_JSON_(dict_, **kwargs)

    @classmethod
    def _from_JSON_(cls, dict_: OrderedDict[str, Any], /, **kwargs):
        return cls(**dict_)


class UserDefined(Identifier):
    def __init__(self, value: Any, class_: str):
        super().__init__(value)
        self.class_ = class_


class Name(Identifier):
    def __init__(self, value: str):
        super().__init__(value)


class SecondaryNames(Identifier):
    def __init__(self, value: list[str]):
        super().__init__(value)


class Abbreviation(Identifier):
    def __init__(self, value: list[str]):
        super().__init__(value)


class SecondaryAbbreviation(Identifier):
    def __init__(self, value: list[str]):
        super().__init__(value)


class CAS(Identifier):
    def __init__(self, value: list[str]):
        super().__init__(value)


class SMILES(Identifier):
    def __init__(self, value: str):
        super().__init__(value)


class ChemicalFormula(Identifier):
    def __init__(self, value: str):
        super().__init__(value)

    def element_count(self, element: str) -> int:
        pattern = re.escape(element) + r'\D*(\d+)'  # look for target, optional non-digits, then digits
        match = re.search(pattern, self.value)
        if match:
            return int(match.group(1))

        # look for target only
        match = re.search(re.escape(element), self.value)
        if match:
            return 1

        return 0


class INCHI(Identifier):
    def __init__(self, value: str):
        super().__init__(value)


class INCHIKey(Identifier):
    def __init__(self, value: str):
        super().__init__(value)
