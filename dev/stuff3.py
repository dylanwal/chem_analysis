from functools import wraps
from typing import Protocol

import numpy as np
from scipy.signal import find_peaks


class SignalProtocol(Protocol):
    x: np.ndarray
    y: np.ndarray


def stuff(a: SignalProtocol):
    return a.x


class A:
    def __init__(self):
        self._x = np.ones(10)
        self._y = np.ones(12)

    @property
    def x(self):
        return np.ones(10)

    @property
    def y(self):
        return np.ones(10)

aa = A()
b = stuff(aa)




import functools

def func(x,y):
    return 1

new_func = functools.partial(func, cen1= 2, cen1=2)

new_func()