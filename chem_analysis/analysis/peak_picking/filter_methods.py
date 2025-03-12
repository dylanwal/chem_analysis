

class HeightFilter(FilterMethods):
    def __init__(self,
                 min_abs: float = None,
                 max_abs: float = None,
                 min_rel: float = None,
                 max_rel: float = None,
                 ):
        """

        Parameters
        ----------
        min_abs
        max_abs
        min_rel
        max_rel
        """
        if min_abs is None and min_abs is None and min_abs is None and max_abs is None:
            raise ValueError("At least on height parameter is required")
        self.min_abs = min_abs
        self.max_abs = max_abs
        self.min_rel = min_rel
        self.max_rel = max_rel

    def run(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        heights = y[index]

        mask = np.ones_like(index, dtype=bool)
        if self.min_abs is not None:
            np.logical_and(heights >= self.min_abs, mask, out=mask)
        if self.max_abs is not None:
            np.logical_and(heights <= self.max_abs, mask, out=mask)
        if self.min_rel is not None:
            np.logical_and(heights >= self.min_rel * np.max(y), mask, out=mask)
        if self.max_rel is not None:
            np.logical_and(heights <= self.max_rel * np.max(y), mask, out=mask)

        return index[mask]


class HeightFilterLocal(FilterMethods):
    def __init__(self,
                 min_abs: float = None,
                 max_abs: float = None,
                 min_rel: float = None,
                 max_rel: float = None,
                 type_: str = 'peak',
                 window: int = 1,
                 ):
        """

        Parameters
        ----------
        min_abs
        max_abs
        min_rel
        max_rel
        type_:
            "peak": will compare to adjacent peaks
            "y": will compare to adjacent y values
        window: int
            will look left and right window size
        """
        if min_abs is None and min_abs is None and min_abs is None and max_abs is None:
            raise ValueError("At least on height parameter is required")
        if not(type_ == 'peak' or type_ == 'y'):
            raise ValueError("Type must be 'peak' or 'y'")
        self.min_abs = min_abs
        self.max_abs = max_abs
        self.min_rel = min_rel
        self.max_rel = max_rel
        self.type_ = type_
        self.window = window

    def run(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        heights = y[index]

        mask = np.ones_like(index, dtype=bool)
        for i in range(index):
            if self.type_ == 'peak':
                window = [max(0, i - self.window), min(len(heights), i + self.window)]
                heights_ = heights[window[0]:window[1]]
            if self.type_ == 'y':
                window = [max(0, i - self.window), min(len(y), i + self.window)]
                heights_ = y[window[0]:window[1]]

            if self.min_abs is not None:
                np.logical_and(heights_ - heights[i] <= self.min_abs, mask, out=mask)
            if self.max_abs is not None:
                np.logical_and(heights_ - heights[i] >= self.min_abs, mask, out=mask)
            if self.min_rel is not None:
                np.logical_and(heights >= self.min_rel * np.max(y), mask, out=mask)
            if self.max_rel is not None:
                np.logical_and(heights <= self.max_rel * np.max(y), mask, out=mask)

        return index[mask]