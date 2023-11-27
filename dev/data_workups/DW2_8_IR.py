
import numpy as np
import chem_analysis as ca


def main():
    data = ca.ir.IRSignalArray.from_file(
        r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-8\DW2_8_IIR2.feather",
    )

    signal = data.get_signal(0)
    fig = ca.ir.plot_signal(signal)
    fig.write_html("temp.html", auto_open=True)


if __name__ == "__main__":
    main()
