import pathlib

import plotly.graph_objs as go

import chem_analysis as ca


def main():
    cal_RI = ca.sec.ConventionalCalibration(lambda time: 10 ** (-0.6 * time + 10.644),
                                    mw_bounds=(160, 1_090_000), name="RI calibration")

    # loading data
    file_path = pathlib.Path(r"data//SEC.csv")
    sec_signal = ca.sec.SECSignal.from_csv(file_path)
    sec_signal.calibration = cal_RI

    # processing
    sec_signal.processor.add(ca.processing.baseline.ImprovedAsymmetricLeastSquared(
        lambda_=1e6,
        p=0.15,
    ))

    # analysis
    peak = ca.analysis.peak_picking.find_peak_largest(sec_signal,
                                                      mask=ca.processing.weigths.Spans((10, 12.2), invert=True)
                                                      )
    result = ca.analysis.integration.rolling_ball(peak, n=45, min_height=0.02, n_points_with_pos_slope=1)

    # plotting
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(sec_signal, fig=fig)
    ca.plotting.calibration(cal_RI, fig=fig)
    ca.plotting.peaks(result, fig=fig)
    fig.data[0].line.color = "blue"
    fig.data[4].fillcolor = "gray"
    # fig.show()
    fig.write_image('figs/sec_data_analysis.png', width=fig.layout.template.layout.width, height=fig.layout.template.layout.height)

    # print results
    print(result.stats_table().to_str())

    # print select results
    table = result.stats_table()
    print("Mn (g/mol):", table.get("mw_n")[0])
    print("D:", table.get("mw_d")[0])


if __name__ == '__main__':
    main()
