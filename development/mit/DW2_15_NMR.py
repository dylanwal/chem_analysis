import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca
from chem_analysis.utils.math import get_slice
from chem_analysis.analysis.line_fitting import fitting_adaptive, DistributionNormal, DistributionMultinomial
from chem_analysis.utils.math import normalize_by_max


def conversion(mca_result) -> tuple[np.ndarray, np.ndarray]:
    return mca_result.C[:, 1] / (mca_result.C[:, 0] + mca_result.C[:, 1]), mca_result.C[:, 3] / (
            mca_result.C[:, 2] + mca_result.C[:, 3])


def plot_mca_results(mcrar, x, D, times, convMA, convDMA):
    times = times - times[0]

    # resulting spectra
    fig1 = go.Figure()
    for i in range(mcrar.ST.shape[0]):
        fig1.add_trace(go.Scatter(y=normalize_by_max(mcrar.ST[i, :]), x=x, name=f"mca_{i}"))
    fig1.add_trace(go.Scatter(y=normalize_by_max(D[0, :]), x=x, name="early"))
    fig1.add_trace(go.Scatter(y=normalize_by_max(D[int(D.shape[0] / 2), :]), x=x, name="middle"))
    fig1.add_trace(go.Scatter(y=normalize_by_max(D[-1, :]), x=x, name="late"))
    fig1.update_layout(autosize=False, width=1200, height=600, font=dict(family="Arial", size=18, color="black"),
                       plot_bgcolor="white", showlegend=True)
    fig1.update_xaxes(title="<b>wavenumber (cm-1) (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray", autorange="reversed")
    fig1.update_yaxes(title="<b>normalized absorbance</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray")

    fig5 = go.Figure()
    for i in range(mcrar.ST.shape[0]):
        fig5.add_trace(go.Scatter(y=normalize_by_max(mcrar.ST[i, :]), x=x, name=f"mca_{i}"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.MA.y_normalized_by_max()), x=pure.MA.x, name="MA"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.PMA.y_normalized_by_max()), x=pure.PMA.x, name="PMA"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.DMAA.y_normalized_by_max()), x=pure.DMAA.x, name="DMAA"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.PDMAA.y_normalized_by_max()), x=pure.PDMAA.x, name="PDMAA"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.DMSO.y_normalized_by_max()), x=pure.DMSO.x, name="DMSO"))
    # fig5.add_trace(go.Scatter(y=normalize_by_max(pure.FL.y_normalized_by_max()), x=pure.FL.x, name="FL"))
    fig5.update_layout(autosize=False, width=1200, height=600, font=dict(family="Arial", size=18, color="black"),
                       plot_bgcolor="white", showlegend=True)
    fig5.update_xaxes(title="<b>wavenumber (cm-1) (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray", autorange="reversed")
    fig5.update_yaxes(title="<b>normalized absorbance</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray")

    fig2 = go.Figure()
    for i in range(mcrar.C.shape[1]):
        fig2.add_trace(go.Scatter(x=times, y=mcrar.C[:, i], name=f"mca_{i}"))
    fig2.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                       plot_bgcolor="white", showlegend=True)
    fig2.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray")
    fig2.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray", range=[0, 1])

    fig3 = go.Figure()
    fig3.add_scatter(x=times, y=convMA, name="convMA")
    fig3.add_scatter(x=times, y=convDMA, name="convDMA")
    fig3.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                       plot_bgcolor="white", showlegend=False)
    fig3.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray")
    fig3.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray", range=[0, 1])

    ca.plotting.plotly_utils.merge_figures([fig1, fig2, fig3, fig5], auto_open=True)


def mca_4(x: np.ndarray, t: np.ndarray, D: np.ndarray):
    C = np.ones((D.shape[0], 4))
    C[:, 0] = 0.5
    C[:, 1] = 0.5
    C[:, 2] = 0.02
    C[:, 3] = 0.02
    C[34, :] = np.array([0.01, 0.01, 0.3, 0.6])
    C[49, :] = np.array([0.20, 0.70, 0.01, 0.01])
    # C[101, :] = np.array([0.20, 0.70, 0.01, 0.01])
    # C[114, :] = np.array([0.01, 0.01, 0.7, 0.29])

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[
            ca.analysis.mca.ConstraintRange(
                [
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                ]
            ),
            ca.analysis.mca.ConstraintConv(index=[0, 1, 2, 3]),
        ],
        st_constraints=[ca.a.mca.ConstraintNonneg()],
        tolerance_increase=10
    )
    mca_result = mca.fit(D, C=C, verbose=True)

    convMA, convDMA = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, convMA, convDMA)
    return np.column_stack((t, convMA, convDMA))


def fit_MA(nmr_array: ca.nmr.NMRSignal2D):
    areas = np.zeros((nmr_array.number_of_signals, 2))
    for i, signal in enumerate(nmr_array.signal_iter()):
        range_ = [3.53, 3.7]
        slice_ = get_slice(signal.x, start=range_[0], end=range_[1])
        x = signal.x[slice_]
        y = signal.y[slice_]
        if np.max(y) < 0.1:
            continue

        result = fitting_adaptive(
            DistributionMultinomial((DistributionNormal, DistributionNormal)),
            x=x, y=y,
            num_trials=20
        )

        # peaks = [
        #     DistributionNormal(x, 1, 3.57, 0.01, scale_bounds= [.01, 1], mean_bounds=[3.54, 3.6], sigma_bounds=[0.002, 0.04]),
        #     DistributionNormal(x, 1, 3.67, 0.01, scale_bounds= [.01, 1], mean_bounds=[3.65, 3.72], sigma_bounds=[0.002, 0.04]),
        # ]
        # result = peak_deconvolution(peaks=peaks, xdata=x, ydata=y)

        PMA = result.multipeak.peaks[0]
        MA = result.multipeak.peaks[1]
        areas[i] = [MA.area(), PMA.area()]

    return areas


def conv_from_normal(nmr_array: ca.nmr.NMRSignal2D):
    MA_areas = fit_MA(nmr_array)

    fig = go.Figure()
    fig.add_scatter(x=nmr_array.y-nmr_array.y[0], y=MA_areas[:, 0], name="MA")
    fig.add_scatter(x=nmr_array.y-nmr_array.y[0], y=MA_areas[:, 1], name="PMA")
    fig.add_scatter(x=nmr_array.y-nmr_array.y[0], y=MA_areas[:, 1]/np.sum(MA_areas, axis=1), name="conv")
    fig.show()

    for i in range(nmr_array.number_of_signals):
        print(nmr_array.y[i], MA_areas[i, 0], MA_areas[i, 1])


def plot_fit(x, y, peaks):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y))
    for peak in peaks:
        fig.add_trace(go.Scatter(x=peak.x, y=peak.y))

    # fig.show()
    fig.write_html("spectra.html", auto_open=True)


def do_integrate(nmr_array: ca.nmr.NMRSignal2D):
    x_ranges = {
        'I_benzene': (7.15, 7.45),  # 6
        'I_MA': (3.6, 3.75),  # 3PMA
        'I_PMA': (3.42, 3.6),  # 3PMA
        'I_DMA_PDMA': (2.60, 3.2),  # 3DMA, 3PDMA
        'I_DMA1': (5.4, 5.75),  # 1 DMA
        'I_DMA2': (6.47, 7.03),  # 1DMA
        'I_DMA_MA': (5.75, 6.4)  # 3 MA, 1 DMA
    }
    results = {key: np.zeros(nmr_array.number_of_signals) for key in x_ranges}
    for i, signal in enumerate(nmr_array.signal_iter()):
        for key, x_range in x_ranges.items():
            results[key][i] = ca.a.integrate.integrate_trapz(signal, x_range=x_range)
    return results


def integrate_array(nmr_array: ca.nmr.NMRSignal2D):
    results = do_integrate(nmr_array)
    I_benzene = results['I_benzene']
    I_MA = results['I_MA']/results['I_benzene']/3
    I_PMA = results['I_PMA']/results['I_benzene']/3
    I_DMA = results['I_DMA1']/results['I_benzene']
    I_DMA2 = results['I_DMA2']/results['I_benzene']
    I_PDMA = (results['I_DMA_PDMA'] - 6*results['I_DMA1'])/results['I_benzene']/6

    x = np.arange(len(I_MA))
    fig = go.Figure()
    smoother = ca.processing.smoothing.Gaussian(sigma=1)
    # _, I_MA = smoother.run_xy(x, I_MA)
    # _, I_PMA = smoother.run_xy(x, I_PMA)
    # _, I_DMA = smoother.run_xy(x, I_DMA)
    # _, I_PDMA = smoother.run_xy(x, I_PDMA)

    fig.add_scatter(x=x, y=I_MA, name="MA")
    fig.add_scatter(x=x, y=I_DMA, name="DMA")
    fig.add_scatter(x=x, y=I_PMA, name="PMA")
    fig.add_scatter(x=x, y=I_DMA2, name="DMA2")
    fig.add_scatter(x=x, y=I_PDMA, name="PDMA")
    fig.show()


    # print
    for i in range(len(I_MA)):
        print(nmr_array.y[i], I_MA[i], I_DMA[i], I_PDMA[i], I_PMA[i])

    # conv_DMA = 1 - I_DMA / (I_DMA + I_PDMA)
    # conv_MA = 1 - I_MA / (I_MA + I_PMA)
    # fig = go.Figure()
    # fig.add_scatter(x=nmr_array.y-1703178980, y=conv_MA)
    # fig.add_scatter(x=nmr_array.y-1703178980, y=conv_DMA)
    # fig.layout.yaxis.range = [0, 1]
    # fig.show()


def main():
    path = pathlib.Path(r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2_15_NMR.feather")
    data = ca.nmr.NMRSignal2D.from_feather(path)
    data.delete([0, 1, 2, 109, 107, 85])

    # process data
    signals = list(data.signal_iter())
    for signal in signals:
        signal.processor.add(
            # ca.processing.translations.AlignMaxValue(range_=(2.5, 2.7), x_value=2.53, temporal_processing=True),
            ca.p.smoothing.Wavelet(threshold=0.015),
            ca.processing.smoothing.Gaussian(sigma=5),
            ca.processing.translations.AlignMaxValue(range_=(7.2, 7.5), x_value=7.33, min_height=0.01),
        )
        print(signal)

    data = ca.nmr.NMRSignal2D.from_signals(signals, unify_method=ca.base.unify_methods.UnifyMethodExpandInterpolate())

    # t_slice = slice(1, None)
    # mask = ca.p.weights.Spans(
    #         [
    #             [2.7, 3.3],
    #             [3.42, 3.8],
    #             [5.4, 7]
    #         ],
    #     )
    # mask = mask.get_mask(data.x, data.y)
    # mca_result_1 = mca_4(data.x[mask], data.y[t_slice], data.z[t_slice, mask])
    # np.savetxt("mca_result.csv", mca_result_1, delimiter=',')

    integrate_array(data)
    # conv_from_normal(data)

    # visualize
    n = 56
    fig = ca.plot.signal(data.get_signal(n, processed=False), raw=True)
    fig = ca.plot.signal(data.get_signal(n), fig=fig)
    # fig = ca.plot.signal2D_slices(data, slice(0,None,5))
    fig.show()

    print("done")


if __name__ == "__main__":
    main()
