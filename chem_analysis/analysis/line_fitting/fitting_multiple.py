

class ResultMultiFitting:
    def __init__(self, results: list[ResultFitting] = None):
        self.results = results or []
        self._best_index: int | None = None

    def __str__(self):
        return f"Best result: {self.best_result}"

    def _determine_best(self):
        if len(self.results) == 0:
            return

        best_score = np.inf
        for i, result in enumerate(self.results):
            if result.best_score < best_score:
                self._best_index = i
                best_score = result.best_score

    @property
    def best_params(self) -> tuple[Any]:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index].best_params

    @property
    def best_score(self) -> PeakModel:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index].model(*self.best_params)

    @property
    def best_result(self) -> ResultFitting:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index]

    @property
    def best_model(self) -> PeakModel:
        if self._best_index is None:
            self._determine_best()
        return self.results[self._best_index].model

    def add_result(self, result: ResultFitting):
        self.results.append(result)
        self._best_index = None



def fitting_multiple_model_adaptive(
        models: Sequence[Type[PeakModel]],
        x: np.ndarray,
        y: np.ndarray,
        trials: Sequence[np.ndarray] = None,
        criteria: Criteria = None,
        full_output: bool = False,
) -> PeakModel | ResultMultiFitting:
    criteria = criteria or BIC

    results = ResultMultiFitting()
    for i, model in enumerate(models):
        fit_result = fitting_adaptive(model, x, y, trials, criteria, full_output=True)
        results.add_result(fit_result)

    check_fit(y, results.best_model(x))
    if full_output:
        return results
    return results.best_model
