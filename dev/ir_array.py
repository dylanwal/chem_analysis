import pathlib

import numpy as np

import chem_analysis as ca


def main():
    file_path = pathlib.Path(r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-5\DW2-5-1-ATIR.csv")
    data = np.loadtxt(file_path, delimiter=",")
    wavenumber = np.flip(data[0, 1:])
    times = data[1:5, 0]
    data = data[1:5, 1:]
    array_ = ca.ir.IRArray(x=wavenumber, y=times, z=data)

    fig = ca.ir.plot_signal(array_.signals[1], x_range=slice(1000, 2000))
    fig.write_html("temp.html", auto_open=True)
    # fig = ca.ir.plot_signal_array_overlap(array_, y_range_index=slice(0, -1, 25))
    # fig.write_html("temp.html", auto_open=True)
    # fig2 = ca.ir.plot_3D_traces(array_, x_range=slice(900, 2000), y_range_index=slice(0, -1, 25))
    # fig2.write_html("temp.html", auto_open=True)

    # mca = ca.algorithms.analysis.multi_component_analysis.MCA()
    # D = array_.z
    # C = np.ones((D.shape[0], 3)) * .5
    # C[0, :] = np.array([1, 0, 0])
    # results = mca.fit(D, C, verbose=True)
    # print(results.C)


if __name__ == "__main__":
    main()
