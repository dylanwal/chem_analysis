import pathlib

import chem_analysis as ca


def main():
    # loading data
    file_path = pathlib.Path(r"data//ir.feather")
    ir_signal = ca.ir.IRSignal.from_feather(file_path)

    # plotting
    fig = ca.plotting.signal(ir_signal)
    fig.show()


if __name__ == '__main__':
    main()
