import os
import pathlib

import chem_analysis as ca


def load_many(path: pathlib.Path):
    files = tuple(os.scandir(path))
    files = sorted(files, key=lambda x: int(x.name[1:]))
    nmr_signals = []
    for file in files:
        nmr_ = ca.nmr.NMRSignal.from_spinsolve_csv(file.path)
        nmr_signals.append(nmr_)
        if len(nmr_signals) == 10:
            break

    return ca.nmr.NMRSignalArray.from_signals(nmr_signals)


def main():
    path = pathlib.Path(r"C:\Users\nicep\Desktop\DW2-7")
    data = load_many(path)
    data.processor.add(ca.processing.re_sampling.CutSpans(x_spans=(-2, 12)))
    print(data)
    # data.raw_data = np.flip(data.raw_data, axis=1)

    signal = data.get_signal(0)
    fig = ca.ir.plot_signal(signal)
    fig.write_html("temp.html", auto_open=True)

    # mca_result_1 = mca_2_mask(data)


if __name__ == "__main__":
    main()
