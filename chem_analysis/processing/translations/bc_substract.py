

class Subtract(BaselineCorrection):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 multiplier: float = 1,
                 ):
        super().__init__(None)
        self.y_sub = y
        self.x_sub = x
        self.multiplier = multiplier

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if len(self.y_sub) == len(y):
            return self.multiplier * self.y_sub

        raise NotImplementedError()  # TODO: x-interpolation

    def get_baseline2D(self, x: np.ndarray, _: np.ndarray, z: np.ndarray) -> np.ndarray:
        baseline = np.empty_like(z)

        for i in range(z.shape[0]):
            baseline[i, :] = self.get_baseline(x, z[i, :])

        return baseline



class SubtractOptimize(BaselineCorrection):
    def __init__(self,
                 y: np.ndarray,
                 x: np.ndarray = None,
                 weights: DataWeight | Iterable[DataWeight] = None,
                 bounds: tuple[float, float] = (-2, 2)
                 ):
        super().__init__(weights)
        self.y_sub = y
        self.x_sub = x
        self.bounds = bounds

        self.multiplier = None

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.weights is not None:
            mask = self.weights.get_mask(x, y)
            x_ = x[mask]
            y_ = y[mask]
            y_sub = self.y_sub[mask]
            x_sub = self.x_sub[mask]
        else:
            x_ = x
            y_ = y
            y_sub = self.y_sub
            x_sub = self.x_sub

        self.multiplier = self._get_multiplier(x_, y_, x_sub, y_sub)
        return self.multiplier * self.y_sub

    def _get_multiplier(self, x: np.ndarray, y: np.ndarray, x_sub: np.ndarray, y_sub: np.ndarray) -> float:
        if len(self.y_sub) == len(y):
            def func(m) -> float:
                return float(np.sum(np.abs(y-m*y_sub)))

            result = minimize_scalar(func, bounds=self.bounds)
            if not result.success:
                raise ValueError(f"'{type(self).__name__}' has not converged.")
            return result.x

        raise NotImplementedError()  # TODO: x-interpolation

    def get_baseline2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        self.multiplier = np.empty_like(y)
        baseline = np.empty_like(z)

        for i in range(z.shape[0]):
            y = z[i, :]
            if self.weights is not None:
                mask = self.weights.get_mask(x, y)
                x_ = x[mask]
                y_ = y[mask]
                y_sub = self.y_sub[mask]
                x_sub = self.x_sub[mask]
            else:
                x_ = x
                y_ = y
                y_sub = self.y_sub
                x_sub = self.x_sub

            multiplier = self._get_multiplier(x_, y_, x_sub, y_sub)
            self.multiplier[i] = multiplier
            baseline[i, :] = multiplier * self.y_sub

        return baseline