
import numpy as np

import chem_analysis as ca

# cal_RI = ca.sec.ConventionalCalibration(lambda time: 10 ** (-0.0107*time**2 + 0.236*time + 5.7072),
#                                         mw_bounds=(589, 2_110_000), name="RI calibration")
cal_RI = ca.sec.ConventionalCalibration(lambda time: 10 ** (-0.2986*time+12.245),
                                        mw_bounds=(589, 2_110_000), name="RI calibration")


def process_one(sig: ca.sec.SECSignal, *, output: bool = False, save_img: bool = False):
    # sig = ca.sec.SECSignal(x_raw=time_, y_raw=y, calibration=cal_RI)
    # sig.processor.add(
    #     ca.p.baseline.BaselineWithMask(
    #         ca.p.baseline.ImprovedAsymmetricLeastSquared(lambda_=1e6, p=0.15),
    #         mask=ca.processing.weigths.Spans(((6, 9), (16.8, 17.2)), invert=True),
    #         save_result=True
    #     )
    # )
    sig.processor.add(
        ca.p.smoothing.Wavelet(threshold=0.5),
        ca.p.smoothing.Gaussian(sigma=60),
        ca.p.baseline.BaselineWithMask(
            ca.p.baseline.Polynomial(degree=0),
            ca.p.weigths.Spans(((20,25)))  # ((22.3,23.4),(30.5,30.6)))
        )
    )

    peaks = ca.analysis.peak_picking.find_peak_largest(sig, mask=ca.processing.weigths.Spans((23, 30)), min_height=10)
    peaks = chem_analysis.analysis.peaks.integration.rolling_ball(peaks, n=150, min_height=0.08, n_points_with_pos_slope=1)
    peaks.peaks = [peak for peak in peaks if peak.properties.area() > 1]

    if output:
        stats_table = peaks.stats_table()
        print(stats_table.write_csv_str(limit_to=["max_x", "mw_d", "mw_n"]))

        fig_base = ca.plotting.signal(sig, raw=True)
        fig_base = ca.plot.signal(sig, fig=fig_base)
        # fig_base = ca.plot.baseline(sig, fig=fig_base)
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
    npzfile = np.load(r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2-15-SEC-DMF-RI.npz")
    x = npzfile['x']
    y = npzfile['y'].T
    array = ca.sec.SECSignalArray(x=x, y=np.arange(y.shape[0]), z=y, calibration=cal_RI)

    # process_one(array.get_signal(64), output=True)
    process_many(list(array.signal_iter()), save_img=True)


if __name__ == "__main__":
    main()

