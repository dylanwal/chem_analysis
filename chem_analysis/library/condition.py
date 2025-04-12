import abc
from typing import Any
from collections import OrderedDict
import inspect

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = inspect.getmembers(class_, lambda a: not (inspect.isroutine(a)))
    list_of_attrs = [a for a in attrs if not (a[0].startswith('_'))]
    return {k: v for k, v in list_of_attrs}


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

    def __str__(self):
        text = str(self.value)
        if self.unit:
            text += f" {self.unit}"
        if self.uncertainty:
            text += f" ({self.uncertainty})"
        return text

    def __repr__(self):
        return self.__str__()

    def to_dict(self, remove_nones: bool = False) -> OrderedDict[str, Any]:
        dict_ = OrderedDict()
        dict_["type"] = type(self).__name__

        attrs = get_class_instance_attributes(self)
        dict_["value"] = attrs.pop('value')
        attrs.pop('unit')

        if remove_nones:
            if self.unit is not None:
                dict_["unit"] = self.unit
            attrs = {k: v for k, v in attrs.items() if v is not None}
        else:
            dict_["unit"] = self.unit

        for k, v in attrs.items():
            dict_[k] = v
        return dict_

    def write_json(self, /, **kwargs) -> OrderedDict[str, Any]:
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


class UserDefined(Condition):
    def __init__(self, value: Any, class_: str):
        super().__init__(value)
        self.class_ = class_


class Method(Condition):
    def __init__(self, value: str):
        super().__init__(value)


class TimeISO(Condition):
    def __init__(self, value: str):
        super().__init__(value)


class Duration(Condition):
    def __init__(self, value: int | float, unit: str = "s", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Temperature(Condition):
    def __init__(self, value: int | float, unit: str = "degC", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Pressure(Condition):
    def __init__(self, value: int | float, unit: str = "kPa", uncertainty: Any | None = None):
        super().__init__(value, unit, uncertainty)


class Solvent(Condition):
    def __init__(self, value: str):
        super().__init__(value)
