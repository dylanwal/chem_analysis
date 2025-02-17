import pathlib
import datetime

import numpy as np
import plotly.graph_objs as go
from scipy.ndimage import gaussian_filter

import chem_analysis as ca

start = datetime.datetime.fromisoformat('2023-12-21T11:16:20')


def main():
    remove = [
        [-594, 0],
        [6576, 7169],
        [13755, 14341],
        [17356, 17726],
        [19995, 20243],
        [21497, 22094],
        [28926, 29167],
    ]

    path = pathlib.Path(r"C:\Users\nicep\Desktop\Book1.csv")
    data = np.loadtxt(path, delimiter=",")
    print(data.shape)

    for r in reversed(remove):  # start from back to front since we are shifting times
        slice_ = ca.utils.math.get_slice(data[:, 0], *r)
        data = np.delete(data, slice_, axis=0)
        data[slice_.start:, 0] -= (r[1] - r[0])  # shift time so gap doesn't exist

    print(data.shape)
    path = path.with_stem(path.stem + "_proc")
    np.savetxt(path, data, delimiter=",")


if __name__ == "__main__":
    main()
