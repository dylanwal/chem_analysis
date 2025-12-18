import pathlib

import plotly.graph_objs as go
import numpy as np
import chem_analysis as ca
from chem_analysis.utils.math import normalize_by_max


class Pure:
    def __init__(self):
        path = r"C:\Users\nicep\Desktop\dynamic_poly\data\IR_standards"
        path = pathlib.Path(path)
        self.MA = ca.ir.IRSignal.from_csv(path / "MA_DMSO_20_percent.csv")
        self.PMA = ca.ir.IRSignal.from_csv(path / "PMA_DMSO_20_percent.csv")
        self.DMA = ca.ir.IRSignal.from_csv(path / "DMAA_DMSO_20_percent.csv")
        self.PDMA = ca.ir.IRSignal.from_csv(path / "PDMAA_DMSO_20_percent.csv")
        self.DMSO = ca.ir.IRSignal.from_csv(path / "DMSO.csv")
        self.FL = ca.ir.IRSignal.from_csv(path / "FL.csv")
        self.mask = None

    def apply_mask(self, mask: np.ndarray):
        self.mask = mask

    @property
    def MA_y(self) -> np.ndarray:
        return self.MA.y[self.mask]

    @property
    def PMA_y(self) -> np.ndarray:
        return self.PMA.y[self.mask]

    @property
    def DMSO_y(self) -> np.ndarray:
        return self.DMSO.y[self.mask]

    @property
    def FL_y(self) -> np.ndarray:
        return self.FL.y[self.mask]


def conversion(mca_result) -> np.ndarray:
    return mca_result.C[:, 1] / (mca_result.C[:, 0] + mca_result.C[:, 1])


def plot_mca_results(mcrar, x, D, times, conv):
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
    pure = Pure()
    for i in range(mcrar.ST.shape[0]):
        fig5.add_trace(go.Scatter(y=normalize_by_max(mcrar.ST[i, :]), x=x, name=f"mca_{i}"))
    fig5.add_trace(go.Scatter(y=normalize_by_max(pure.MA.y_normalized_by_max()), x=pure.MA.x, name="MA"))
    fig5.add_trace(go.Scatter(y=normalize_by_max(pure.PMA.y_normalized_by_max()), x=pure.PMA.x, name="PMA"))
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
    fig3.add_trace(go.Scatter(x=times, y=conv))
    fig3.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                       plot_bgcolor="white", showlegend=False)
    fig3.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray")
    fig3.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                      linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                      gridwidth=1, gridcolor="lightgray", range=[0, 1])

    ca.plotting.plotly_utils.merge_figures([fig1, fig2, fig3, fig5], auto_open=True)


def mca_2(x: np.ndarray, t: np.ndarray, D: np.ndarray, pure: Pure):
    C = np.ones((D.shape[0], 2)) * .5
    C[0, :] = np.array([.7, 0.3])

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[ca.analysis.mca.ConstraintNonneg(), ca.analysis.mca.ConstraintConv()],
        st_constraints=[],
        tolerance_increase=100
    )
    mca_result = mca.fit(D, C=C, verbose=True)

    conv = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, conv)
    return np.column_stack((t, conv))


def mca_4(x: np.ndarray, t: np.ndarray, D: np.ndarray, pure: Pure):
    # C = mca_4_ST(x, t, D, pure)
    # C = C[t]

    C = np.ones((D.shape[0], 4))
    C[:, 0] = 0.5
    C[:, 1] = 0.5
    C[:, 2] = 0.02
    C[:, 3] = 0.02
    C[0, :] = np.array([1, 0.01, 0.01, 0.01])
    C[40, :] = np.array([0.6, 0.4, 0.01, 0.01])

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[
            ca.analysis.mca.ConstraintConv(index=[0, 1]),
            ca.analysis.mca.ConstraintRange(
                [
                    (0, 1),
                    (0, 1),
                    (-0.2, 1),
                    (0, 3),
                ]
            )
        ],
        st_constraints=[],
        tolerance_increase=100
    )
    mca_result = mca.fit(D, C=C, verbose=True)

    conv = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, conv)
    return np.column_stack((t, conv))


def mca_2_ST(x: np.ndarray, t: np.ndarray, D: np.ndarray, pure: Pure):
    ST = np.ones((2, D.shape[1]))
    ST[0, :] = pure.MA
    ST[1, :] = pure.PMA

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[ca.analysis.mca.ConstraintNonneg(), ca.analysis.mca.ConstraintConv()],
        st_constraints=[],
        tolerance_increase=100
    )
    mca_result = mca.fit(D, ST=ST, st_fix=[0,1], verbose=True)

    conv = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, conv)
    return np.column_stack((t, conv))


def mca_5_ST(x: np.ndarray, t: np.ndarray, D: np.ndarray, pure: Pure):
    ST = np.ones((5, D.shape[1]))
    ST[0, :] = pure.MA_y
    ST[1, :] = pure.PMA_y
    ST[2, :] = pure.DMSO_y * 0.2
    ST[3, :] = pure.FL_y
    noise = D[-1]
    noise = noise / np.max(noise)
    mask = noise > 0.2
    noise[mask] = .2
    mask = noise < -0.2
    noise[mask] = -0.2
    ST[4, :] = noise

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[
            ca.analysis.mca.ConstraintConv(index=[0, 1]),
            ca.analysis.mca.ConstraintRange(
                [
                    (0, 1),
                    (0, 1),
                    (-3, 0),
                    (0, 3),
                    (0.01, 0.3)
                ]
            )
        ],
        st_constraints=[ca.analysis.mca.ConstraintNonneg(index=[4])],
        tolerance_increase=100
    )
    mca_result = mca.fit(D, ST=ST, st_fix=[0, 1, 2, 3, 4], verbose=True)

    conv = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, conv)
    return np.column_stack((t, conv))


def mca_4_ST(x: np.ndarray, t: np.ndarray, D: np.ndarray, pure: Pure):

    ST = np.ones((4, D.shape[1]))
    ST[0, :] = pure.MA_y
    ST[1, :] = pure.PMA_y
    ST[2, :] = pure.DMSO_y * 0.2
    ST[3, :] = pure.FL_y

    print("working")
    mca = ca.analysis.mca.MultiComponentAnalysis(
        max_iters=200,
        c_constraints=[
            ca.analysis.mca.ConstraintConv(index=[0, 1]),
            ca.analysis.mca.ConstraintRange(
                [
                    (0, 1),
                    (0, 1),
                    (-0.2, 1),
                    (0, 3),
                ]
            )
        ],
        st_constraints=[],
        tolerance_increase=100
    )
    mca_result = mca.fit(D, ST=ST, st_fix=[0, 1, 2, 3], verbose=True)

    conv = conversion(mca_result)
    plot_mca_results(mca_result, x, D, t, conv)
    return np.column_stack((t, conv))


def main():
    data = ca.ir.IRSignal2D.from_feather(
        r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2_15_IR.feather"
    )
    # data.processor.add(ca.processing.baselines.SubtractOptimize(data.data[-2, :]))
    # data.to_feather(r"C:\Users\nicep\Desktop\post_doc_2022\Data\polymerizations\DW2-14\DW2_14_IR_proc.feather")

    # signal = data.get_signal(300)
    # fig = ca.plot.signal(signal)
    # fig.add_trace(go.Scatter(x=signal.x_raw, y=signal.y_raw))
    # fig.write_html("temp.html", auto_open=True)

    # data.processor.add(ca.processing.re_sampling.CutOffValue(x_span=529, cut_off_value=0.015))
    # data.processor.add(ca.processing.translations.ScaleMax(range_=(1700, 1800)))
    # data.processor.add(ca.processing.smoothing.GaussianTime(sigma=3))
    # data.processor.add(ca.processing.smoothing.ExponentialTime(a=0.7))
    # data.processor.add(ca.processing.baselines.Polynomial(
    #         degree=1,
    #         weights=ca.processing.weigths.AdaptiveDistanceMedian(threshold=0.2)
    # ))
    # data.processor.add(
    #     ca.processing.baselines.Polynomial(
    #         degree=1,
    #         weights=ca.processing.weigths.Spans(x_spans=(1900, 2000), invert=True)
    #     )
    # )
    # data.processor.add(ca.processing.smoothing.Gaussian(sigma=2))
    # data.processor.add(ca.processing.translations.ScaleMax(range_=(1700, 1800)))

    ## get conversion
    mca_result_1 = mca_2(*mca_pre(data))
    for i in range(mca_result_1.shape[0]):
        print(mca_result_1[i, 0], mca_result_1[i, 1])

    ## making gifs
    # create_gif(data)
    # create_gif_surface(data)


def mca_pre(data: ca.ir.IRSignal2D):
    pure = Pure()
    t_slice = slice(60, None)
    mask = ca.processing.weigths.Slices(
        [
            ca.utils.math.get_slice(data.x, start=None, end=760),
            ca.utils.math.get_slice(data.x, start=875, end=1100),
            ca.utils.math.get_slice(data.x, start=1350, end=1600),
            ca.utils.math.get_slice(data.x, start=1900, end=None),
        ],
    )
    mask = mask.get_mask(data.x, data.y)
    pure.apply_mask(mask)

    return data.x[mask], data.y[t_slice], data.z[t_slice, mask], pure


if __name__ == "__main__":
    main()





