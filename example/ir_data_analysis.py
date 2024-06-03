import pathlib

import chem_analysis as ca


def main():
    # loading data
    file_path = pathlib.Path(r"data//ir.feather")
    ir_signal = ca.ir.IRSignal.from_feather(file_path)

    # plotting
    fig = ca.plotting.signal(ir_signal)
    # fig.show()
    fig.write_image('figs/ir_data_analysis.png', width=fig.layout.template.layout.width, height=fig.layout.template.layout.height)


if __name__ == '__main__':
    main()
