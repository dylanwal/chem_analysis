import numpy as np


def left_pad(x: np.ndarray, pad_amount: int = 5):
    delta = x[1] - x[0]
    return np.linspace(x[0] - delta * pad_amount, x[0] - delta, pad_amount)


def right_pad(x: np.ndarray, pad_amount: int = 5):
    delta = x[-1] - x[-2]
    return np.linspace(x[-1] + delta, x[-1] + delta * pad_amount, pad_amount)


def get_x_pad(x: np.ndarray, pad_amount: int = 5, side: str = "both"):
    FLAG_FLIP = False
    if x[0] > x[-1]:
        x = np.flip(x)
        FLAG_FLIP = True

    if side == "left" or side == "both":
        left_side = left_pad(x, pad_amount)
    else:
        left_side = []
    if side == "right" or side == "both":
        right_side = right_pad(x, pad_amount)
    else:
        right_side = []

    x_new = np.concatenate((left_side, x, right_side))
    if FLAG_FLIP:
        return np.flip(x_new)
    return x_new


def pad_edges_flat_with_x(x: np.ndarray, y: np.ndarray, pad_amount: int = 5) -> tuple[np.ndarray, np.ndarray]:
    y_new = pad_edges_flat(y, pad_amount)
    x_new = get_x_pad(x, pad_amount)
    return x_new, y_new


def pad_edges_flat(y: np.ndarray, pad_amount: int = 5) -> np.ndarray:
    left_side = [y[0]] * np.ones(pad_amount)
    right_side = [y[-1]] * np.ones(pad_amount)
    return np.concatenate((left_side, y, right_side))


def pad_edges_polynomial_with_x(
        x: np.ndarray,
        y: np.ndarray,
        *,
        pad_amount: int = 5,
        degree: int = 1,
        window: int = 5,
        side: str = "both"
) -> tuple[np.ndarray, np.ndarray]:
    """

    Parameters
    ----------
    y
    x
    pad_amount
    degree
    window
    side:
        "both" = both sides
        "left" = low side (left side)
        "right" = high side (right side)

    Returns
    -------

    """
    FLAG_FLIP = False
    if x[0] > x[-1]:
        x = np.flip(x)
        y = np.flip(y)
        FLAG_FLIP = True

    if side == "left" or side == "both":
        left_pad_, left_side = pad_edges_polynomial_single(y, x=x, pad_amount=pad_amount, degree=degree, window=window)
    else:
        left_side, left_pad_ = [], []
    if side == "right" or side == "both":
        right_pad_, right_side = pad_edges_polynomial_single(y, x=x, pad_amount=pad_amount, degree=degree,
                                                             window=window, left_side=False)
    else:
        right_side, right_pad_ = [], []

    x_new = np.concatenate((left_pad_, x, right_pad_))
    y_new = np.concatenate((left_side, y, right_side))
    if FLAG_FLIP:
        return np.flip(x_new), np.flip(y_new)

    return x_new, y_new


def pad_edges_polynomial(
        y: np.ndarray,
        *,
        pad_amount: int = 5,
        degree: int = 1,
        window: int = 5,
        side: str = "both"
) -> np.ndarray:
    """

    Parameters
    ----------
    y
    pad_amount
    degree
    window
    side:
        "both" = both sides
        "left" = low side (left side)
        "right" = high side (right side)

    Returns
    -------

    """
    if side == "left" or side == "both":
        _, left_side = pad_edges_polynomial_single(y, pad_amount=pad_amount, degree=degree, window=window)
    else:
        left_side = []
    if side == "right" or side == "both":
        _, right_side = pad_edges_polynomial_single(y, pad_amount=pad_amount, degree=degree, window=window,
                                                    left_side=False)
    else:
        right_side = []

    return np.concatenate((left_side, y, right_side))


def pad_edges_polynomial_single(
        y: np.ndarray,
        *,
        x: np.ndarray = None,
        pad_amount: int = 5,
        degree: int = 1,
        window: int = 5,
        left_side: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """

    Parameters
    ----------
    y
    x
    pad_amount
    degree
    window
    left_side:
        0 = low side (left side)
        1 = high side (right side)

    Returns
    -------

    """
    if pad_amount < 1:
        raise ValueError("Pad length needs to be greater than 1")

    if left_side:  # left
        if x is None:
            x = np.arange(len(y))
        x_pad = left_pad(x, pad_amount)
        x_window = x[:window]
        y_window = y[:window]

    else:  # right
        if x is None:
            x = np.arange(len(y))
        x_pad = right_pad(x, pad_amount)
        x_window = x[-window:]
        y_window = y[-window:]

    poly = np.polynomial.Polynomial.fit(x_window, y_window, degree)
    return x_pad, poly(x_pad)


def run_local():
    x = np.linspace(0, 10, 80)
    y = np.sin(x)

    x_new, y_new = pad_edges_polynomial_with_x(x, y)
    print(y_new)

    import plotly.graph_objs as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name="base"))
    fig.add_trace(go.Scatter(x=x_new, y=y_new, mode="lines+markers", name="new"))
    fig.show()


if __name__ == "__main__":
    run_local()
