from __future__ import annotations
from typing import Sequence, TypeVar
import abc

import numpy as np
from scipy.special import voigt_profile

from chem_analysis.utils.math import rescale_array


class PeakModel(abc.ABC):
    def __str__(self):
        args_text = ','.join([arg + f": {getattr(self, arg):0.3f}" for arg in self.__slots__])
        return f"{type(self).__name__}({args_text})"

    @abc.abstractmethod
    def __call__(self, x: np.ndarray) -> np.ndarray:
        ...

    def number_of_params(self) -> int:
        return len(self.__slots__)

    @classmethod
    def initial_guess_generator(cls,
                                x: np.ndarray,
                                y: np.ndarray,
                                trials: Sequence[str] | str | None = None,
                                number_trials: int = 5,
                                ) -> Sequence[np.ndarray]:
        raise NotImplementedError()


class DistributionNormal(PeakModel):
    __slots__ = ("scale", "mean", "sigma")

    def __init__(self,
                 scale: int | float,
                 mean: int | float,
                 sigma: int | float,
                 ):
        super().__init__()
        self.scale = scale
        self.mean = mean
        self.sigma = sigma

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.scale / (self.sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - self.mean) ** 2 / (2 * self.sigma ** 2))

    @classmethod
    def initial_guess_generator(cls,
                                x: np.ndarray,
                                y: np.ndarray,
                                trials: Sequence[str] | str | None = None,
                                num_trials: int = 5,
                                ) -> Sequence[np.ndarray]:
        if trials is None:
            trials = cls.__slots__
        else:
            if isinstance(trials, str):
                trials = [trials]
            if not all([t not in cls.__slots__ for t in trials]):
                raise ValueError(f"Invalid trial provided. "
                                 f"\n\tprovided trials: {trials} "
                                 f"\n\taccepted trials: {cls.__slots__}."
                                 )

        from scipy.stats import qmc
        sampler = qmc.LatinHypercube(d=len(trials))
        sample = sampler.random(num_trials)

        counter = 0
        if "scale" in trials:
            y_max = y.max()
            bounds = (y_max-0.2*(y_max-y.min()), y_max)
            scale = rescale_array(sample[:, counter], bounds[0], bounds[1])
            counter += 1
        else:
            scale = np.ones(num_trials) * np.max(y)

        if "mean" in trials:
            span = x.max() - x.min()
            bounds = (x.min()+0.15*span, x.max()-0.15*span)  # 0.15 is to move bounds more center
            mean = rescale_array(sample[:, counter], bounds[0], bounds[1])
            counter += 1
        else:
            mean = np.ones(num_trials) * np.mean(x)

        if "sigma" in trials:
            span = x.max() - x.min()
            bounds = (span*0.01, span)
            sigma = rescale_array(sample[:, counter], bounds[0], bounds[1])
        else:
            sigma = np.ones(num_trials) * 0.5*(x.max() - x.min())

        return [np.array([scale[i], mean[i], sigma[i]]) for i in range(num_trials)]


class DistributionNormal2(PeakModel):
    __slots__ = ("scale", "mean", "sigma")

    def __init__(self,
                 scale1: int | float,
                 mean1: int | float,
                 sigma: int | float,
                 ):
        super().__init__()
        self.scale = scale
        self.mean = mean
        self.sigma = sigma

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.scale / (self.sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - self.mean) ** 2 / (2 * self.sigma ** 2))


class DistributionCauchy(PeakModel):
    __slots__ = ("scale", "mean", "gamma")

    def __init__(self,
                 scale: int | float,
                 mean: int | float,
                 gamma: int | float,
                 ):
        super().__init__()
        self.scale = scale
        self.mean = mean
        self.gamma = gamma

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.scale / (np.pi * self.gamma * (1 + ((x - self.mean) / 2) ** 2))


class DistributionVoigt(PeakModel):
    __slots__ = ("scale", "mean", "gamma")

    def __init__(self,
                 scale: int | float,
                 mean: int | float,
                 sigma: int | float,
                 gamma: int | float
                 ):
        """
        gamma = 0 normal
        sigma = 0 cauchy distribution
        """
        super().__init__()
        self.scale = scale
        self.mean = mean
        self.sigma = sigma
        self.gamma = gamma

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.scale * voigt_profile(x - self.mean, sigma=self.sigma, gamma=self.gamma)


class DistributionMultinomial(PeakModel):
    def __init__(self, models: list[PeakModel]):
        super().__init__()
        self.models = models

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return np.sum(tuple(model(x) for model in self.models))

    def initial_guess_generator(self,
                                x: np.ndarray,
                                y: np.ndarray,
                                trials: Sequence[str] = None,
                                num_trials: int = 5,
                                ) -> Sequence[np.ndarray]:
        model_trials = [model.initial_guess_generator(x, y, trials, num_trials) for model in self.models]

        arr = np.array(model_trials)
        trials = arr.swapaxes(0,1).reshape(20,6)
        return trials

    def set_params(self, *args, **kwargs):
        return self