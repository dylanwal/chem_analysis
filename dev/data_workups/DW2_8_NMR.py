
import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


# def intgration_ben(signals: list[ca.nmr.NMRSignal]):
#     for sig in signals:
#         peaks = ca.analysis.peak_picking.scipy_find_peaks(sig)
#         if len(peak)
#         result = ca.analysis.boundary_detection.rolling_ball(peaks)
#         result.get_stats()


def main():
    path = pathlib.Path(r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-8\DW2_8_NMR2.feather")
    nmr_array = ca.nmr.NMRSignalArray.from_file(path)
    nmr_array.processor.add(ca.processing.translations.AlignMax(range_=(2.2, 2.7)))

    # I_bond = ca.analysis.integrate.integrate(nmr_array, x_range=(6.11, 6.4))
    # mask = I_bond < 0.05
    # I_benzene = ca.analysis.integrate.integrate(nmr_array, x_range=(7.2, 7.5))
    # I_ratio = I_bond/I_benzene
    # I_ratio[mask] = np.max(I_ratio)
    # I_poly = ca.analysis.integrate.integrate(nmr_array, x_range=(3.52, 3.71))
    # mask = I_poly < 0.05
    # I_poly[mask] = 0
    # I_mon = ca.analysis.integrate.integrate(nmr_array, x_range=(3.71, 3.83))
    # # mask = I_mon < 0.05
    # # I_mon[mask] = 0
    # for i in range(len(nmr_array.time_zeroed)):
    #    print(nmr_array.time_zeroed[i], I_poly[i]/(I_mon[i]+I_poly[i]))
    # fig = go.Figure(go.Scatter(x=nmr_array.time_zeroed, y=1 - I_ratio/np.max(I_ratio)))
    # fig.add_trace(go.Scatter(x=nmr_array.time_zeroed, y=I_poly/(I_mon+I_poly),
    # text=[f'X: {x}, Index: {i}' for i, x in enumerate(nmr_array.time_zeroed)],
    # hoverinfo='text',
    #                          ))
    # fig.write_html("conv.html", auto_open=True)

    signal = nmr_array.get_signal(5)
    fig = ca.ir.plot_signal(signal)
    fig.write_html("temp.html", auto_open=True)



    #signals = [nmr_array.get_signal(i) for i in range(len(nmr_array.time))]


if __name__ == "__main__":
    main()
