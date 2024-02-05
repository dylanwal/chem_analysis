import pathlib

import chem_analysis as ca


def load_pure():
    path = r"C:\Users\nicep\Desktop\post_doc_2022\Data\polymerizations"
    path = pathlib.Path(path)
    MA = ca.ir.IRSignal.from_file(path / "ATIR_MA.csv")
    PMA = ca.ir.IRSignal.from_file(path / "ATIR_PMA.csv")
    DMSO = ca.ir.IRSignal.from_file(path / "ATIR_DMSO.csv")
    FL = ca.ir.IRSignal.from_file(path / "ATIR_perflourohexane.csv")

    return MA, PMA, DMSO, FL


def main():
    MA, PMA, DMSO, FL = load_pure()

    PMA.processor.add(
        ca.processing.baselines.AdaptiveAsymmetricLeastSquared(
            lambda_=1e7))

    fig = ca.plotting.signal(PMA)
    fig = ca.plotting.signal_raw(PMA, fig=fig)
    fig = ca.plotting.baseline(PMA, fig=fig)
    fig.show()

    for i in range(len(PMA.x)):
        print(f"{PMA.x[i]}, {PMA.y[i]}")


if __name__ == "__main__":
    main()
