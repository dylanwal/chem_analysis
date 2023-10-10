from __future__ import annotations
import abc

import numpy as np


def get_subclasses(cls, depth: int | None = 0, _count: int = 0) -> set[type]:
    subs = set(cls.__subclasses__())
    if _count < depth or depth is None:
        return subs.union(*(get_subclasses(i, depth, _count) for i in subs))

    return subs


class ProcessingMethod(abc.ABC):

    @classmethod
    def processing_classes(cls) -> set[type]:
        return get_subclasses(cls, depth=0)

    @classmethod
    def processing_algorithms(cls) -> set[type]:
        return get_subclasses(cls, depth=None)

    @abc.abstractmethod
    def run(self, *args, **kwargs):
        ...


class Processor(list):
    """
    Processor
    """
    def __repr__(self):
        return f"Processor: {len(self)} methods"

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        for method in self:
            x, y = method.run(x, y)
        return x, y
