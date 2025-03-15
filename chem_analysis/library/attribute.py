import abc
import base64
from typing import Any
from collections import OrderedDict

import numpy as np

from chem_analysis.library.condition import Condition
from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def numpy_to_JSON(array: np.ndarray, encoding: str = "list") -> str:
    if encoding == "binary":
        return f"{array.shape}|{array.dtype}|" + base64.b64encode(array.tobytes()).decode('ASCII')
    if encoding == "list":
        return array.tolist()

    raise ValueError("encoding not supported.")


def parse_shape(shape: str) -> list[int]:
    return [int(i) for i in shape[1:-1].split(",") if i]


def JSON_to_numpy(input_: str, encoding: str = "list") -> np.ndarray:
    if encoding == "binary":
        shape,  dtype, data = input_.split("|", maxsplit=2)
        return np.frombuffer(base64.b64decode(input_), dtype=dtype).reshape(parse_shape(shape))
    if encoding == "list":
        return np.array(input_)

    raise ValueError("numpy encoding not supported.")


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = class_.__dict__.items()
    return {k: v for k, v in list(attrs) if v is not None}


class Attribute(MixinSubClassList, abc.ABC):
    __slots__ = "value", "unit", "uncertainty", "conditions"
    MINI_KEY: str

    def __init__(self,
                 value: Any,
                 unit: str | None = None,
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        self.value = value
        self.unit = unit
        self.uncertainty = uncertainty
        if isinstance(conditions, Condition):
            conditions = [conditions]
        self.conditions = conditions or []

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
        return class_._from_JSON_(dict_, *args, **kwargs)

    @classmethod
    def _from_JSON_(cls, dict_: OrderedDict[str, Any], *args, **kwargs):
        return cls(**dict_)


class MolarMass(Attribute):
    MINI_KEY = "mw"

    def __init__(self,
                 value: float | int,
                 unit: str = 'g/mol',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Color(Attribute):
    MINI_KEY = "color"

    def __init__(self,
                 value: str,
                 conditions: Condition = None
                 ):
        super().__init__(value, conditions=conditions)


class BoilingTemperature(Attribute):
    MINI_KEY = "btemp"

    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class MeltingTemperature(Attribute):
    MINI_KEY = "mtemp"

    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Density(Attribute):
    MINI_KEY = "den"

    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class MassSpectrum(Attribute):
    MINI_KEY = "ms"

    def __init__(self,
                 value: np.ndarray,  # [m,2]  first col is m/z; second is intensity.
                 unit: str = 'Da',
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, conditions=conditions)

    def to_json(self, numpy_encoding: str = "list") -> OrderedDict[str, Any]:
        dict_ = self.to_dict()
        dict_["value"] = numpy_to_JSON(dict_["value"])
        return dict_

    @classmethod
    def _from_JSON_(cls, dict_: OrderedDict[str, Any], **kwargs):
        dict_["value"] = JSON_to_numpy(dict_["value"], encoding=kwargs["numpy_encoding"])
        return cls(**dict_)


class RetentionTime(Attribute):
    MINI_KEY = "rt"

    def __init__(self,
                 value: float | int,
                 unit: str = 'min',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class ResponseFactor(Attribute):
    MINI_KEY = "rt"

    def __init__(self,
                 value: float | int,
                 unit: str = 'min',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class VaporPressure(Attribute):
    MINI_KEY = "vp"

    def __init__(self,
                 value: float | int,
                 unit: str = 'kPa',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Solubility(Attribute):
    MINI_KEY = "vp"

    def __init__(self,
                 value: float | int,
                 unit: str = 'g/ml',
                 uncertainty: Any | None = None,
                 conditions: Condition = None
                 ):
        super().__init__(value, unit, uncertainty, conditions)
