

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

cal_RI = ca.sec.ConventionalCalibration(lambda time: 10 ** (-0.6 * time + 10.644),
                                        mw_bounds=(160, 1_090_000), name="RI calibration")


def process_one(sig: ca.sec.SECSignal, *, output: bool = False, save_img: bool = False):
    # sig = ca.sec.SECSignal(x_raw=time_, y_raw=y, calibration=cal_RI)
    sig.processor.add(
        ca.p.baseline.BaselineWithMask(
            ca.p.baseline.ImprovedAsymmetricLeastSquared(lambda_=1e6, p=0.15),
            mask=ca.processing.weigths.Spans(((6, 9), (16.8, 17.2)), invert=True),
            save_result=True
        )
    )

    peaks = ca.analysis.peak_picking.find_peak_largest(sig, mask=ca.processing.weigths.Spans((10, 12.2), invert=True))
    peaks = ca.analysis.integration.rolling_ball(peaks, n=45, min_height=0.05, n_points_with_pos_slope=1)
    peaks.peaks = [peak for peak in peaks if peak.area() > 1]

    if output:
        stats_table = peaks.stats_table()
        print(stats_table.write_csv_str())

        fig_base = ca.plotting.signal(sig, raw=True)
        fig_base = ca.plot.signal(sig, fig=fig_base)
        fig_base = ca.plot.baseline(sig, fig=fig_base)
        fig_base.layout.yaxis.range = (-1, 50)

        fig = ca.plot.signal(sig)
        fig = ca.plot.peaks(peaks, fig=fig)
        fig = ca.plot.calibration(sig.calibration, fig=fig)
        fig.layout.yaxis.range = (-1, 50)

        ca.plotting.plotly_utils.merge_figures([fig_base, fig], auto_open=True)

    if save_img:
        fig = ca.plot.signal(sig)
        fig = ca.plot.peaks(peaks, fig=fig)
        fig = ca.plot.calibration(sig.calibration, fig=fig)
        fig.layout.yaxis.range = (-1, 50)
        fig.write_image(f"img/signal{sig.id_}.png")

    return peaks.stats_table()


def process_many(signals: list[ca.sec.SECSignal], *, output: bool = False, save_img: bool = False):
    table = None
    for sig in signals:
        if table is None:
            table = process_one(sig, output=output, save_img=save_img)
        else:
            table.join(process_one(sig, output=output, save_img=save_img))
        print(f"{sig.id_} done")

    print(table.write_csv_str(limit_to=["max_x", "mw_d", "mw_n"]))


def main():
    npzfile = np.load(r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2-15-SEC-RI.npz")
    x = npzfile['x']
    y = npzfile['y'].T
    array = ca.sec.SECSignalArray(x=x, y=np.arange(y.shape[0]), z=y, calibration=cal_RI)

    # process_one(array.get_signal(10), output=True)
    process_many(list(array.signal_iter()), save_img=True)


if __name__ == "__main__":
    main()

