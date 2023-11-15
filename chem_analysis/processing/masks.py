import abc
from typing import Iterable

import numpy as np


class DataMask(abc.ABC):

    @abc.abstractmethod
    def get_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class DataMaskSlice(DataMask):
    def __init__(self, slices: slice | Iterable[slice]):
        self.slices = slices

    def get_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(self.slices, list):
            data_x = []
            for slice_ in self.slices:
                data_x.append(x[slice_])
            data_y = []
            for slice_ in self.slices:
                data_y.append(y[slice_])
            return np.vstack(data_x), np.vstack(data_y)

        return x[self.slices], y[self.slices]


class DataMaskMultiPoint(DataMask):
    def __init__(self, indexes: Iterable):
        self.indexes = indexes

    def get_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x[self.indexes], y[self.indexes]


def get_index_distance_from_median(y: np.ndarray, m: float = 2):
    dist_from_median = np.abs(y - np.median(y))
    median_deviation = np.median(dist_from_median)
    if median_deviation != 0:
        scale_distances_from_median = dist_from_median / median_deviation
        return scale_distances_from_median > m  # True is an outlier

    return np.zeros_like(y)  # no outliers


class DataMaskDistanceMedian(DataMask):
    def __init__(self, norm_distance: float = 2):
        self.norm_distance = norm_distance
        self.indexes = None

    def get_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self.indexes = get_index_distance_from_median(y, self.norm_distance)
        return x[self.indexes], y[self.indexes]


def get_distance_remove(y: np.ndarray, amount: float = 0.7, speed: float = 0.1, max_rounds: int = 1000):
    stop_len = len(y) * amount
    indexes = np.ones_like(y, dtype=bool)
    for i in range(max_rounds):
        dist_from_median = np.abs(y - np.median(y[indexes]))
        median_deviation = np.median(dist_from_median[indexes])
        if median_deviation == 0:
            median_deviation = np.mean(dist_from_median[indexes])
        scale_distances_from_median = dist_from_median / median_deviation
        keep_points = scale_distances_from_median < np.max(scale_distances_from_median[indexes]) * (1-0.5*np.exp(-i/speed))
        if np.sum(keep_points) == 0:
            break
        indexes = np.bitwise_and(indexes, keep_points)

        if np.sum(indexes) < stop_len:
            break

    return indexes


class DataMaskDistanceRemove(DataMask):
    def __init__(self, amount: float = 0.5, speed: float = 10, max_rounds: int = 1000, keep_ends: bool = True):
        self.amount = amount
        self.speed = speed
        self.max_rounds = max_rounds
        self.keep_ends = keep_ends
        self.indexes = None

    def get_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self.indexes = get_distance_remove(y, self.amount, self.speed, self.max_rounds)
        if self.keep_ends:
            self.indexes[0] = True
            self.indexes[-1] = True
        return x[self.indexes], y[self.indexes]
