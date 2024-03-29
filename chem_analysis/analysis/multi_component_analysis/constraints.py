import abc
from typing import Iterable

import numpy as np


class Constraint(abc.ABC):
    """ Abstract class for constraints

    Parameters
    ----------
    copy : bool
        Make copy of input data, A; otherwise, overwrite (if mutable)
    """

    def __init__(self, copy=True):
        self.copy = copy

    @abc.abstractmethod
    def transform(self, A):
        """ Transform A input based on constraint """


class ConstraintNonneg(Constraint):
    """
    Non-negativity constraint. All negative entries made 0.
    """

    def __init__(self, index=None, copy=False):
        """ A must be non-negative"""
        super().__init__(copy)
        self.index = index

    def transform(self, A):
        """ Apply nonnegative constraint"""
        if self.index is not None:
            if self.copy:
                A = np.copy(A)
            selected_columns = A[:, self.index]
            A[:, self.index] = selected_columns * (selected_columns > 0)
            return A

        if self.copy:
            return A * (A > 0)
        else:
            A *= (A > 0)
            return A


class ConstraintCumsumNonneg(Constraint):
    """
    Cumulative-Summation non-negativity constraint. All negative
     entries made 0.

    """

    def __init__(self, axis=-1, copy=False):
        """ A must be non-negative"""
        super().__init__(copy)
        self.axis = axis

    def transform(self, A):
        """ Apply cumsum nonnegative constraint"""

        if self.copy:
            return A * (np.cumsum(A, self.axis) > 0)
        else:
            A *= (np.cumsum(A, self.axis) > 0)
            return A


class ConstraintZeroEndPoints(Constraint):
    """
    Enforce the endpoints (or the mean over a range) is zero

    """

    def __init__(self, axis=-1, span=1, copy=False):
        """

        Parameters
        ----------
        axis : int
            Axis to operate on

        span : int
            Number of pixels along the ends to average.
        copy
        """
        super().__init__(copy)
        if [0, 1, -1].count(axis) != 1:
            raise TypeError('Axis must be 0, 1, or -1')

        self.axis = axis
        self.span = span

    def transform(self, A):
        """ Apply cumsum nonnegative constraint"""

        pix_vec = np.arange(A.shape[self.axis])
        if (self.axis == 0):
            if self.span == 1:
                slope = (A[-1, :] - A[0, :]) / (pix_vec[-1] - pix_vec[0])
                intercept = A[0, :]
            else:
                slope = ((A[-self.span:, :].mean(axis=0) -
                          A[:self.span, :].mean(axis=0)) /
                         (pix_vec[-self.span:].mean() -
                          pix_vec[:self.span].mean()))
                intercept = (A[:self.span, :] -
                             np.dot(pix_vec[:self.span, None],
                                     slope[None, :])).mean(axis=0)
            if self.copy:
                return A - np.dot(pix_vec[:, None],
                                   slope[None, :]) - intercept[None, :]
            else:
                A -= (np.dot(pix_vec[:, None],
                              slope[None, :]) + intercept[None, :])
                return A
        else:
            if self.span == 1:
                slope = (A[:, -1] - A[:, 0]) / (pix_vec[-1] - pix_vec[0])
                intercept = A[:, 0]
            else:
                slope = ((A[:, -self.span:].mean(axis=1) -
                          A[:, :self.span].mean(axis=1)) /
                         (pix_vec[-self.span:].mean() -
                          pix_vec[:self.span].mean()))
                intercept = (A[:, :self.span] -
                             np.dot(slope[:, None],
                                     pix_vec[None, :self.span])).mean(axis=1)

            if self.copy:
                return A - np.dot(slope[:, None],
                                   pix_vec[None, :]) - intercept[:, None]
            else:
                A -= (np.dot(slope[:, None],
                              pix_vec[None, :]) + intercept[:, None])
                return A


class ConstraintZeroCumSumEndPoints(Constraint):
    """
    Enforce the endpoints of the cumsum (or the mean over a range) is
    near-zero. Note: this is an approximation.


    """

    def __init__(self, nodes=None, axis=-1, copy=False):
        """

        Parameters
        ----------
        copy : bool
            Make copy of input data, A; otherwise, overwrite (if mutable)
        nodes : list of int
            In addition to end-points, other points to ensure are approximately 0
        axis : int
            Axis to operate on

        span : int
            Number of pixels along the ends to average.
        """
        super().__init__(copy)

        self.nodes = nodes
        if [0, 1, -1].count(axis) != 1:
            raise TypeError('Axis must be 0, 1, or -1')

        self.axis = axis

    def transform(self, A):
        """ Apply cumsum nonnegative constraint"""

        meaner = A.mean(self.axis)

        if self.nodes:
            self.nodes = set(self.nodes)
            self.nodes.update({0, A.shape[self.axis]})
            self.nodes.discard(-1)
            self.nodes = list(self.nodes)

            if (self.axis == 0):
                if self.copy:
                    temp = 1.0 * A
                    for num in range(len(self.nodes) - 1):
                        n0 = self.nodes[num]
                        n1 = self.nodes[num + 1]
                        temp[n0:n1, :] -= temp[n0:n1, :].mean(self.axis)[None, :]
                    return temp
                else:
                    for num in range(len(self.nodes) - 1):
                        n0 = self.nodes[num]
                        n1 = self.nodes[num + 1]
                        A[n0:n1, :] -= A[n0:n1, :].mean(self.axis)[None, :]
                    return A
            else:
                if self.copy:
                    temp = 1.0 * A
                    for num in range(len(self.nodes) - 1):
                        n0 = self.nodes[num]
                        n1 = self.nodes[num + 1]
                        temp[:, n0:n1] -= temp[:, n0:n1].mean(self.axis)[:, None]
                    return temp
                else:
                    for num in range(len(self.nodes) - 1):
                        n0 = self.nodes[num]
                        n1 = self.nodes[num + 1]
                        A[:, n0:n1] -= A[:, n0:n1].mean(self.axis)[:, None]
                    return A
        else:
            if (self.axis == 0):
                if self.copy:
                    return A - meaner[None, :]
                else:
                    A -= meaner[None, :]
                    return A
            else:
                if self.copy:
                    return A - meaner[:, None]
                else:
                    A -= meaner[:, None]
                    return A


class ConstraintNorm(Constraint):
    """
    Normalization constraint.

    """

    def __init__(self, axis=-1, fix=None, copy=False):
        """

        Parameters
        ----------
        axis : int
            Which axis of input matrix A to apply normalization acorss.
        fix : list
            Keep fix-axes as-is and normalize the remaining axes based on the
            residual of the fixed axes.
        set_zeros_to_feature : int
            Set all samples which sum-to-zero across axis to 1 for a particular
             feature (See Notes)
        copy : bool
            Make copy of input data, A; otherwise, overwrite (if mutable)

        Notes
        -----

        -   For set_zeros_to_feature, assuming the data represents concentration
             with a matrix [n_samples, n_features] and the axis is across the
             features, for every sample that sums to 0 across axis, would be
             replaced with a vector [n_features] of zeros except at
             set_zeros_to_feature, which would equal 1. I.e., this pixel is
             now pure substance of index value set_zeros_to_feature.

        """
        super().__init__(copy)
        if fix is None:
            self.fix = fix
        elif isinstance(fix, int):
            self.fix = [fix]
        elif isinstance(fix, (list, tuple)):
            self.fix = fix
        elif isinstance(fix, np.ndarray):
            if np.issubdtype(fix.dtype, np.integer):
                self.fix = fix.tolist()
            else:
                raise TypeError('ndarrays must be of dtype int')
        else:
            raise TypeError('Parameter fix must be of type None, int, list,',
                            'tuple, ndarray')

        if not ((axis == 0) | (axis == 1) | (axis == -1)):
            raise ValueError('Axis must be 0,1, or -1')
        self.axis = axis

    def transform(self, A):
        """ Apply normalization constraint """

        if self.copy:
            if self.axis == 0:
                if not self.fix:  # No fixed axes
                    return A / A.sum(axis=self.axis)[None, :]
                else:  # Fixed axes
                    fix_locs = self.fix
                    not_fix_locs = [v for v in np.arange(A.shape[0]).tolist()
                                    if self.fix.count(v) == 0]
                    scaler = np.ones(A.shape)
                    div = A[not_fix_locs, :].sum(axis=0)[None, :]
                    div[div == 0] = 1
                    scaler[not_fix_locs, :] = ((1 - A[fix_locs, :].sum(axis=0)[None, :]) / div)

                    return A * scaler
            else:  # Axis = 1 / -1
                if not self.fix:  # No fixed axes
                    return A / A.sum(axis=self.axis)[:, None]
                else:  # Fixed axis
                    fix_locs = self.fix
                    not_fix_locs = [v for v in np.arange(A.shape[-1]).tolist()
                                    if self.fix.count(v) == 0]
                    scaler = np.ones(A.shape)
                    div = A[:, not_fix_locs].sum(axis=-1)[:, None]
                    div[div == 0] = 1
                    scaler[:, not_fix_locs] = ((1 - A[:, fix_locs].sum(axis=-1)[:, None]) / div)

                    return A * scaler
        else:  # Overwrite original data
            if A.dtype != float:
                raise TypeError('A.dtype must be float for',
                                'in-place math (copy=False)')

            if self.axis == 0:
                if not self.fix:  # No fixed axes
                    A /= A.sum(axis=self.axis)[None, :]
                    return A
                else:  # Fixed axes
                    fix_locs = self.fix
                    not_fix_locs = [v for v in np.arange(A.shape[0]).tolist()
                                    if self.fix.count(v) == 0]
                    scaler = np.ones(A.shape)
                    div = A[not_fix_locs, :].sum(axis=0)[None, :]
                    div[div == 0] = 1
                    scaler[not_fix_locs, :] = ((1 - A[fix_locs, :].sum(axis=0)[None, :]) / div)
                    A *= scaler
                    return A
            else:  # Axis = 1 / -1
                if not self.fix:  # No fixed axes
                    A /= A.sum(axis=self.axis)[:, None]
                    return A
                else:  # Fixed axis
                    fix_locs = self.fix
                    not_fix_locs = [v for v in np.arange(A.shape[-1]).tolist()
                                    if self.fix.count(v) == 0]
                    scaler = np.ones(A.shape)
                    div = A[:, not_fix_locs].sum(axis=-1)[:, None]
                    div[div == 0] = 1
                    scaler[:, not_fix_locs] = ((1 - A[:, fix_locs].sum(axis=-1)[:, None]) / div)
                    A *= scaler
                    return A


class ConstraintReplaceZeros(Constraint):
    """
    Samples that sum-to-zero across axis are replaced with a vector of 0's except
    for a 1 at feature if a single value. In a concentration context, e.g., samples with
    no concentration are replaced with 100% concentration of a set feature. If multiple
    features given, equal amounts of each feature (summing to 1) are used.



    """

    def __init__(self, axis=-1, feature=None, fval=1, copy=False):
        """
            Parameters
        ----------
        axis : int
            Which axis of input matrix A to apply normalization acorss.
        feature : int, list, tuple
            Set all samples which sum-to-zero across axis to fval for a particular
             feature (or fractional) for multiple features.
        fval : float
            Value of summation across axis of replacement vector.
        copy : bool
            Make copy of input data, A; otherwise, overwrite (if mutable)

        """
        super().__init__(copy)
        self.fval = fval
        if feature is None:
            self.feature = feature
        elif isinstance(feature, int):
            self.feature = [feature]
        elif isinstance(feature, (list, tuple)):
            self.feature = feature
        elif isinstance(feature, np.ndarray):
            if np.issubdtype(feature.dtype, np.integer):
                self.feature = feature.tolist()
            else:
                raise TypeError('ndarrays must be of dtype int')
        else:
            raise TypeError('Parameter feature must be of type None, int, list, tuple, ndarray')

        if not ((axis == 0) | (axis == 1) | (axis == -1)):
            raise ValueError('Axis must be 0,1, or -1')
        self.axis = axis

    def transform(self, A):
        """ Apply constraint """

        if self.feature:
            replacement = np.zeros(A.shape[self.axis])
            replacement[self.feature] = self.fval
            replacement /= replacement.sum()
            replacement *= self.fval

            if self.copy:
                A_out = 1 * A
                if self.axis == 0:
                    A_out[:, A_out.sum(axis=0) == 0] = replacement[:, None]
                else:  # Axis 1 / -1
                    A_out[A_out.sum(axis=-1) == 0] = replacement
                return A_out
            else:
                if self.axis == 0:
                    A[:, A.sum(axis=0) == 0] = replacement[:, None]
                else:  # Axis 1 / -1
                    A[A.sum(axis=-1) == 0] = replacement
                return A
        else:
            return A


class ConstraintPlanarize(Constraint):
    """
    Set a particular target to a plane


    """

    def __init__(self, target, shape, use_vals_above=None, use_vals_below=None,
                 lims_to_plane=True, scaler=None, recalc_scaler=False, copy=False):
        """

        Parameters
        ----------

        target : int, list, tuple
            Target numbers to set to a fitted plane
        shape : tuple, list
            Shape of array (M,N) which is (Y,X)
        use_vals_above : float
            Only calculate based on values above (not including)
        use_vals_below : float
            Only calculate based on values below (not including)
        lims_to_plane : bool
            The returned plane will be limited to the range of the optionally supplied
            use_vals_below, use_vals above.
        scaler : float
            A large value that is much bigger than any values in the input array.
            Needed to ensure SVD properly creates plane. If None, auto-calculates.
        recalc_scaler : bool
            Auto-calculate for every new input (does not use previously provided or
            calculated value)
        copy : bool
            Make copy of input data, A; otherwise, overwrite (if mutable)

        Notes
        -----

        -   This uses an SVD to calculate the vector normal to the plane that fits the input data. It
            assumes that the 3rd singular vector is the normal; thus, the x and y vectors for the data need
            be larger than the variance of the input data. Scaler enables this by scaling the auto-generated
            x and y vectors to be much larger than the max-min of the input data
        """
        super().__init__(copy)
        if isinstance(target, int):
            self.target = [target]
        elif isinstance(target, (list, tuple, np.ndarray)):
            self.target = target
        else:
            raise TypeError('target must be an int, list, 2D ndarray, or tuple')

        self.shape = shape
        self.copy = copy
        self.scaler = scaler
        self.recalc = recalc_scaler
        self.use_above = use_vals_above
        self.use_below = use_vals_below
        self.lims_to_plane = lims_to_plane

        self._x = None
        self._y = None
        self._X = None
        self._Y = None

        if scaler is not None:
            self._setup_xy(scaler)

    def _setup_xy(self, scaler):

        self.scaler = scaler
        self._x = scaler * np.arange(self.shape[1], dtype=float)
        self._y = scaler * np.arange(self.shape[0], dtype=float)

        self._X, self._Y = np.meshgrid(self._x, self._y)
        self._X = self._X.ravel()
        self._Y = self._Y.ravel()

    def transform(self, A):
        """ Set targets, t, to fit planes """

        if (self.scaler is None) | (self.recalc):
            self._setup_xy(1e3 * np.abs(A.max() - A.min()))

        if self.copy:
            A2 = 1 * A

        for t in self.target:
            X2 = 1 * self._X
            Y2 = 1 * self._Y
            Z2 = 1 * A[:, t]

            Stack = np.vstack((X2, Y2, Z2))

            if self.use_above is not None:
                X2 = X2[Z2 > self.use_above]
                Y2 = Y2[Z2 > self.use_above]
                Z2 = Z2[Z2 > self.use_above]
                Stack = np.vstack((X2, Y2, Z2))

            if self.use_below is not None:
                X2 = X2[Z2 < self.use_below]
                Y2 = Y2[Z2 < self.use_below]
                Z2 = Z2[Z2 < self.use_below]
                Stack = np.vstack((X2, Y2, Z2))

            Stack = Stack - Stack.mean(axis=-1)[:, None]

            U, s, Vh = np.linalg.svd(Stack, full_matrices=False)
            norm_to_plane = 1 * U[:, -1]

            plane = (((-norm_to_plane[0] * (self._X - X2.mean())) -
                      (norm_to_plane[1] * (self._Y - Y2.mean()))) /
                     norm_to_plane[2]) + Z2.mean()

            if self.lims_to_plane:
                if self.use_above is not None:
                    plane[plane < self.use_above] = self.use_above
                if self.use_below is not None:
                    plane[plane > self.use_below] = self.use_below
            if not self.copy:
                A[:, t] = plane
            else:
                A2[:, t] = plane
        if self.copy:
            return A2
        else:
            return A


class ConstraintConv(Constraint):
    """
    Conversion constraint. sum(C) = 1
    """

    def __init__(self, index=None, copy=False):
        super().__init__(copy)
        self.index = index

    def transform(self, A):
        """ Apply nonnegative constraint"""
        if self.index is not None:
            if self.copy:
                A = np.copy(A)
            selected_columns = A[:, self.index]
            row_sums = selected_columns.sum(axis=1)
            selected_columns /= row_sums[:, np.newaxis]
            A[:, self.index] = selected_columns
            return A

        row_sums = A.sum(axis=1)
        if self.copy:
            return A / row_sums[:, np.newaxis]
        else:
            A /= row_sums[:, np.newaxis]
            return A


class ConstraintRange(Constraint):
    """
    values must stay within range
    """

    def __init__(self, range_: Iterable[tuple[float | None, float | None] | None] = None, copy=False):
        """

        Parameters
        ----------
        range_:
            None are skipped
            example: [
                None,
                None,
                [0, 1]
                [None, 1]
                [1, None]
            ]
        copy
        """
        super().__init__(copy)
        self.range_ = range_

    def transform(self, A):
        """ Apply nonnegative constraint"""
        if self.copy:
            A = np.copy(A)

        for i, r in enumerate(self.range_):
            if r is not None:
                if r[0] is not None:
                    mask = A[:, i] < r[0]
                    A[mask, i] = r[0]
                if r[1] is not None:
                    mask = A[:, i] > r[1]
                    A[mask, i] = r[1]

        return A
