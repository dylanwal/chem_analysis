from typing import Sequence, Any, Type
import logging
import warnings

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from sklearn.metrics import r2_score

from chem_analysis.analysis.line_fitting.peak_models import PeakModel
from chem_analysis.analysis.line_fitting.fitting_criteria import Criteria, BIC

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=OptimizeWarning)


class ResultFit:
    def __init__(self, params: tuple[Any], score: float):
        self.params = params
        self.score = score

    def __str__(self):
        params = ", ".join(f'{p:0.3f}' for p in self.params)
        return f"Score: {self.score:0.2f} | Params: {params}"

    def __repr__(self):
        return self.__str__()


class ResultFitting:
    def __init__(self, model: Type[PeakModel], results: list[ResultFit] = None):
        self.model = model
        self.results = results or []
        self._best_index: int | None = None

    def __str__(self):
        return f"Model: {self.model.__name__} | Best score: {self.best_score:0.2f} | Best params: {self.best_params}"

    def __repr__(self):
        return self.__str__()
        
    def _determine_best(self):
        if len(self.results) == 0:
            return
        
        best_score = np.inf
        for i, result in enumerate(self.results):
            if result.score < best_score:
                self._best_index = i
                best_score = result.score
    
    @property
    def best_params(self) -> tuple[Any]:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index].params

    @property
    def best_score(self) -> float:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index].score
    
    def add_result(self, result: ResultFit):
        self.results.append(result)
        self._best_index = None


def check_fit(y: np.ndarray, y_new: np.ndarray):
    r2 = r2_score(y, y_new)
    if r2 < 0.95:
        logger.warning('fit was not good enough')
        
    
def fitting_simple(
        model: Type[PeakModel],
        x: np.ndarray, 
        y: np.ndarray, 
        initial_guess: Sequence[Any] = None
):
    kwargs = {}
    if initial_guess is not None:
        kwargs['initial'] = initial_guess
    popt, _ = curve_fit(model, x, y, **kwargs)
    check_fit(y, model(*popt)(x))
    return popt
    

def fitting_adaptive(
        model: Type[PeakModel],
        x: np.ndarray,
        y: np.ndarray,
        trials: Sequence[np.ndarray] | int = None,
        num_trials: int = 5,
        criteria: Criteria = None,
        full_output: bool = False,
) -> tuple | ResultFitting | None:
    trials = trials if isinstance(trials, Sequence) else model.initial_guess_generator(x, y, trials, num_trials)
    criteria = criteria or BIC
        
    results = ResultFitting(model)
    for i, trial in enumerate(trials):
        try:
            def func(x_: np.ndarray, *args):
                return model(*args)(x_)
            popt, _ = curve_fit(func, x, y, p0=trial)
            model_ = model(*popt)
            results.add_result(
                ResultFit(
                    tuple(popt),
                    criteria(x, y, model_(x), model_.number_of_params())
                )
            )
            
        except RuntimeError:
            continue  # Skip failed fits

    check_fit(y, model(*results.best_params)(x))
    if full_output:
        return results
    return results.best_params




