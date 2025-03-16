import numpy as np
from scipy.linalg import lstsq


def deconvolution(
        combined_signal: np.ndarray,
        pure_signals: np.ndarray,
        normalize_signal: bool = True,
        normalize_output: bool = True,
) -> np.ndarray:
    if normalize_signal:
        # normalize to total ion count
        pure_signals = pure_signals.T / np.sum(pure_signals, axis=1)
        combined_signal = combined_signal/np.sum(combined_signal)

    weights, residuals, rank, s = lstsq(pure_signals, combined_signal)

    if normalize_output:
        return weights / np.sum(weights)
    return weights


def local_run():
    combined_signal = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    pure_signals = np.array([
        [1, 2, 3, 4, 5, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 6, 7, 8, 9, 10]
    ])

    weights = deconvolution(combined_signal, pure_signals)
    print("Optimal weights:", weights)


    # combined_signal = np.array([1, 3, 3, 1, 0, 8, 9, 0, 0, 0], dtype=np.float64)
    # pure_signals = np.array([
    #     [1.1, 1, 0, 1.2, 0, 0, 0, 0, 0, 0],
    #     [0, 1, 1, 0, 0, 8, 9, 0, 0, 0]
    # ], dtype=np.float64)
    #
    # weights = deconvolution(combined_signal, pure_signals)
    # print("Optimal weights:", weights)

    pure_signals = (pure_signals.T/np.sum(pure_signals, axis=1)).T
    combined_signal = combined_signal/np.sum(combined_signal)

    import plotly.graph_objects as go
    fig = go.Figure()
    y = np.dot(weights, pure_signals)
    x = list(range(1, 11))
    fig.add_bar(x=x, y=combined_signal)
    fig.add_bar(x=x, y=y)
    fig.add_bar(x=x, y=weights[0]*pure_signals[0, :])
    fig.add_bar(x=x, y=weights[1]*pure_signals[1, :])
    fig.show()


if __name__ == '__main__':
    local_run()
