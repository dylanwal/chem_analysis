import abc
from typing import Any
from collections import OrderedDict


from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = class_.__dict__.items()
    return {k: v for k, v in list(attrs) if v is not None}


class Identifier(MixinSubClassList, abc.ABC):
    __slots__ = "value"
    MINI_KEY: str

    def __init__(self, value: Any):
        self.value = value

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

    def to_dict(self) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["type"] = type(self).__name__
        dict_["value"] = self.value

        attrs = get_class_instance_attributes(self)
        attrs.pop('value')
        for attr in attrs:
            dict_[attr[0]] = attr[1]
        return dict_

    def to_json(self, *args, **kwargs) -> OrderedDict[str, Any]:
        return self.to_dict()

    @classmethod
    def _from_JSON(cls, dict_: OrderedDict[str, Any], *args, **kwargs):
        class_ = dict_.pop("type")
        for k in cls.sub_classes():
            if k.__name__ == class_:
                class_ = k
        return class_._from_JSON_(**dict_)

    @classmethod
    def _from_JSON_(cls, dict_: OrderedDict[str, Any], *args, **kwargs):
        return cls(**dict_)


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


class INCHI(Identifier):
    def __init__(self, value: str):
        super().__init__(value)


class INCHIKey(Identifier):
    def __init__(self, value: str):
        super().__init__(value)
