import logging
import difflib
from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np


@dataclass
class UpToDate(Protocol):
    _up_to_date: bool


class ProcessStep:
    def __init__(self, func: Callable, kwargs: dict = None):
        self.name = func.__name__
        self.func = func
        self.kwargs = kwargs

    def __repr__(self):
        return f"{self.name}; kwargs:{self.kwargs.keys()}"

    def __call__(self, *args, **kwargs):
        if kwargs is not None:
            kwargs = {**kwargs, **self.kwargs}

        return self.func(*args, **kwargs)

    def run(self, *args, **kwargs):
        return self.func(*args, **kwargs, **self.kwargs)


class Pipeline:
    """ Pipline

    Attributes
    ----------
    steps: list[ProcessStep]
        List of ProcessStep
    limit: int
        maximum number of reference allowed
        default is 1000

    Methods
    -------
    add(item)
        Adds item to reference list
    remove(item)
        Removes item from reference list
    """

    def __init__(self, steps: (Callable, list[Callable]) = None, limit: int = 1000, _logger=None,
                 up_to_date: UpToDate = None):
        """
        Parameters
        ----------
        steps:
            objects will be passed to self.add()
        limit: int
            maximum number of reference allowed
            default is 10000
        """
        self._steps = []
        self.limit = limit
        self.count = 0

        if _logger is None:
            self._logger = logging
        else:
            self._logger = _logger

        if steps is not None:
            self.add(steps)

        self.up_to_date = up_to_date

    def __repr__(self):
        if self.count < 4:
            return "; ".join([repr(obj) for obj in self.steps])

        return "; ".join([repr(obj) for obj in self.steps[:2]]) + "; ..."

    def __call__(self):
        return self.steps

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.steps[item]
        elif isinstance(item, str):
            index = self._get_index_from_name(item)
            return self.steps[index]
        elif isinstance(item, slice):
            return [self.steps[i] for i in range(*item.indices(len(self.steps)))]
        else:
            mes = f"{item} not found."
            self._logger.error(mes)
            raise ValueError(mes)

    def __len__(self):
        return len(self._steps)

    def __iter__(self):
        for obj in self._steps:
            yield obj

    @property
    def steps(self):
        return self._steps

    def _get_index_from_name(self, item: str) -> int:
        """ get index from name

        Given an item name, return item index in list.
        This matching is done with difflib, so slight typing errors won't result in errors.

        Parameters
        ----------
        item: str
            name of item you are trying to find

        Returns
        -------
        index: int
            index of item in self._reference list

        Raises
        ------
        Exception
            If item name not found.
        """
        values = [i.name for i in self._steps]
        text = difflib.get_close_matches(word=item, possibilities=values, n=1, cutoff=0.8)
        if text:
            return values.index(text[0])
        else:
            mes = f"'{item}' not found."
            self._logger.error(mes)
            raise ValueError(mes)

    def add(self, steps: (Callable, list[Callable]), **kwargs):
        """ Add
        Adds steps to reference list.
        * if one process step is provided, kwargs are the kwargs for that ProcessStep

        * if you use partial(func) then you can pass a list of functions


        Parameters
        ----------
        steps:
            object that you want to add to the list

        Raises
        -------
        Exception
            If invalid object is provided. An object that does not lead to a valid reference.
        """
        if self.up_to_date is not None:
            self.up_to_date._up_to_date = False

        if isinstance(steps, list):
            if all([isinstance(step, Callable) for step in steps]):
                steps = [ProcessStep(step) for step in steps]
                self._steps += steps
                self.count += len(steps)
                return
            else:
                raise ValueError(f"Invalid list of steps provide: {[str(type(step)) for step in steps]}")

        if isinstance(steps, Callable):
            step = ProcessStep(steps, kwargs)
            self._steps.append(step)
            self.count += 1
            return

        raise TypeError(f"expected: Callable, received: {type(steps)} ")

    def remove(self, steps):
        """ Remove
        Removes object from reference list.

        Parameters
        ----------
        steps:
            object that you want to remove to the list
        """
        if self.up_to_date is not None:
            self.up_to_date._up_to_date = False

        if not isinstance(steps, list):
            steps = [steps]

        remove = []
        for step in steps:
            if isinstance(step, (str, int, slice)):
                step = self[step]

            if steps not in self._steps:
                self._logger.error(f"'{self._steps}'is not in list, so it can't be removed.")
                continue

            remove.append(step)

        if not remove:
            return

        # loop through 'remove list' to remove objs
        for step in remove:
            self._steps.remove(step)
            self.count -= 1

    def as_dict(self) -> list:
        """ Returns list of references for serialization."""
        return [obj.as_dict() for obj in self.steps]

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        for step in self.steps:
            x, y = step.run(x, y)
        return x, y
