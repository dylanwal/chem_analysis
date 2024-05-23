import pytest

import numpy as np

import chem_analysis.utils.math as utils_math


def test_map_discrete_x_axis():
    a = utils_math.map_discrete_x_axis(np.arange(15), np.array([1, 4, 8, 9]), np.array([1, 16, 64, 81]))
    assert np.all(a == np.array([0, 1, 0, 0, 16, 0, 0, 0, 64, 81, 0, 0, 0, 0, 0]))
