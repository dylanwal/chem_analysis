import abc
from typing import Any
from collections import OrderedDict

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = class_.__dict__.items()
    return {k: v for k, v in list(attrs) if v is not None}


class Condition(MixinSubClassList, abc.ABC):
    __slots__ = "value", "unit", "uncertainty"
    MINI_KEY: str

    def __init__(self,
                 value: Any,
                 unit: str | None = None,
                 uncertainty: Any | None = None,
                 ):
        self.value = value
        self.unit = unit
        self.uncertainty = uncertainty

    def to_dict(self) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["type"] = type(self).__name__
        dict_["value"] = self.value
        dict_["unit"] = self.unit

        attrs = get_class_instance_attributes(self)
        attrs.pop('value')
        attrs.pop('unit')
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


class Method(Condition):
    MINI_KEY = "method"

    def __init__(self, value: str):
        super().__init__(value)


class TimeISO(Condition):
    MINI_KEY = "iso"

    def __init__(self, value: str):
        super().__init__(value)


class Duration(Condition):
    MINI_KEY = "dur"

    def __init__(self, value: int | float, unit: str = "s", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Temperature(Condition):
    MINI_KEY = "temp"

    def __init__(self, value: int | float, unit: str = "degC", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Pressure(Condition):
    MINI_KEY = "pres"

    def __init__(self, value: int | float, unit: str = "kPa", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Solvent(Condition):
    MINI_KEY = "sol"

    def __init__(self, value: str):
        super().__init__(value)
