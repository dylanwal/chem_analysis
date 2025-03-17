import abc
import base64
from typing import Any
from collections import OrderedDict
import inspect

import numpy as np

from chem_analysis.library.condition import Condition
from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def numpy_to_JSON(array: np.ndarray, encoding: str = "list") -> str:
    if encoding == "binary":
        return f"b'{','.join(str(i) for i in array.shape)}|{array.dtype}|" + base64.b64encode(array.tobytes()).decode('ASCII')
    if encoding == "list":
        return array.tolist()

    raise ValueError("encoding not supported.")


def JSON_to_numpy(input_: str, encoding: str = "list") -> np.ndarray:
    if encoding == "binary" or input_.startswith("b'"):
        input_ = input_.replace("b'", "")
        shape,  dtype, data = input_.split("|", maxsplit=2)
        shape = [int(i) for i in shape.split(",")]
        return np.frombuffer(base64.b64decode(data), dtype=dtype).reshape(shape)
    if encoding == "list":
        return np.array(input_)

    raise ValueError("numpy encoding not supported.")


def get_class_instance_attributes(class_) -> dict[str, Any]:
    attrs = inspect.getmembers(class_, lambda a: not (inspect.isroutine(a)))
    list_of_attrs = [a for a in attrs if not (a[0].startswith('_'))]
    return {k: v for k, v in list_of_attrs}


class Attribute(MixinSubClassList, abc.ABC):
    __slots__ = "value", "unit", "uncertainty", "conditions"
    MINI_KEY: str

    def __init__(self,
                 value: Any,
                 unit: str | None = None,
                 uncertainty: int | float | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        self.value = value
        self.unit = unit
        self.uncertainty = uncertainty
        if isinstance(conditions, Condition):
            conditions = [conditions]
        self.conditions = conditions or []

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
            attrs = {k: v for k, v in attrs.items() if v}
        else:
            dict_["unit"] = self.unit

        for k, v in attrs.items():
            if isinstance(v, (list, tuple)) and len(v) >= 1:
                v_ = [0]*len(v)
                for i, value in enumerate(v):
                    if hasattr(value, "to_dict"):
                        v_[i] = value.to_dict(remove_nones=remove_nones)
                    else:
                        v_[i] = value
                v = v_
            dict_[k] = v
        return dict_

    def to_json(self, /, **kwargs) -> OrderedDict[str, Any]:
        dict_ = self.to_dict(remove_nones=True)
        if isinstance(self.value, np.ndarray):
            dict_["value"] = numpy_to_JSON(dict_["value"], kwargs.get("numpy_encoding", "list"))
        return dict_

    @classmethod
    def _from_JSON(cls, dict_: OrderedDict[str, Any], /, **kwargs):
        class_ = dict_.pop("type")
        for k in cls.sub_classes():
            if k.__name__ == class_:
                class_ = k
        if "conditions" in dict_:
            dict_["conditions"] = [Condition._from_JSON(cond, **kwargs) for cond in dict_["conditions"]]
        return class_._from_JSON_(dict_, **kwargs)

    @classmethod
    def _from_JSON_(cls, dict_: OrderedDict[str, Any], /, **kwargs):
        if isinstance(dict_["value"], str) and dict_["value"].startswith("b'"):
            dict_["value"] = JSON_to_numpy(dict_["value"], encoding=kwargs.get("numpy_encoding", "binary"))
        if isinstance(dict_["value"], list):
            try:
                dict_["value"] = np.array(dict_["value"])
            except ValueError:
                pass

        return cls(**dict_)


class UserDefined(Attribute):
    def __init__(self,
                 value: Any,
                 class_: str,
                 unit: str = None,
                 uncertainty: int | float | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)
        self.class_ = class_


class MolarMass(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'g/mol',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Color(Attribute):
    def __init__(self,
                 value: str,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, conditions=conditions)


class BoilingTemperature(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class MeltingTemperature(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Density(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'degC',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class MassSpectrum(Attribute):
    def __init__(self,
                 value: np.ndarray,  # [m,2]  first col is m/z; second is intensity.
                 unit: str = 'Da',
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, conditions=conditions)


class RetentionTime(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'min',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class ResponseFactor(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'min',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class VaporPressure(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'kPa',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)


class Solubility(Attribute):
    def __init__(self,
                 value: float | int,
                 unit: str = 'g/ml',
                 uncertainty: Any | None = None,
                 conditions: Condition | list[Condition] = None,
                 ):
        super().__init__(value, unit, uncertainty, conditions)
