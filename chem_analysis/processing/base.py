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


class Processor:
    """
    Processor
    """
    def __init__(self, methods: list[ProcessingMethod] = None):
        self._methods: list[ProcessingMethod] = [] if methods is None else methods
        self.processed = False

    def __repr__(self):
        return f"Processor: {len(self)} methods"

    def __len__(self):
        return len(self._methods)

    @property
    def methods(self) -> list[ProcessingMethod]:
        return self._methods

    def add(self, *args: ProcessingMethod):
        self._methods += args
        self.processed = False

    def insert(self, index: int, method: ProcessingMethod):
        self._methods.insert(index, method)
        self.processed = False

    def delete(self, method: int | ProcessingMethod):
        if isinstance(method, ProcessingMethod):
            self._methods.remove(method)
        else:
            self._methods.pop(method)
        self.processed = False

    def run(self, x: np.ndarray, y: np.ndarray, z: np.ndarray | None = None) \
            -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
        for method in self._methods:
            if z is None:
                x, y = method.run(x, y)
            else:
                x, y, z = method.run(x, y, z)

        self.processed = True
        if z is None:
            return x, y
        return x, y, z
