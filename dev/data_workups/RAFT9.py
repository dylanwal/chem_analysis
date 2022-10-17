import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc
import chem_analysis.analysis.utils.plot_format as chem_pf


def main():
    # calibration
    cal_RI = chem.Calibration(lambda time: 10 ** (-0.6 * time + 10.644), lb=160, ub=1_090_000, name="RI calibration")

    # load data and create signals
    file_name = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Instrument\polymerization\RAFT9-L20-F100uL_min.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    # df = df.drop(["2per", "5per", "10per_2", "20per_2", "30per_2", "60per.1", "40per.1"], axis=1)
    signals = []
    for col in df.columns:
        signals.append(chem.SECSignal(name=col, ser=df[col], cal=cal_RI,
                                      x_label="retention time (min)", y_label="signal"))

    # process signals
    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7, deg=2)
        sig.peak_picking(lb=10.7, ub=11.8)

    # plotting
    colors = chem_pf.get_plot_color(len(signals))
    fig = signals[-1].plot(auto_open=False, op_cal=False, normalize=True, color=colors[0])
    op_cal = False
    for sig, color in zip(reversed(signals[:-1]), colors[1:]):
        if color == colors[-1]:
            op_cal = True
        fig = sig.plot(fig, auto_open=False, op_cal=op_cal, normalize=True, color=color)
    # fig.update_layout(legend=dict(x=.5, y=.95))
    # fig.data[4].fillcolor = 'rgba(172,24,25,0.5)'
    fig.update_xaxes(range=[8, 16])
    fig.write_html('temp.html', auto_open=True)

    # print out table
    df_table = signals[0].stats_df()
    for sig in signals[1:]:
        df_table = pd.concat([df_table, sig.stats_df()], axis=1)

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    print(df_table)


def main2():
    # calibration
    cal_RI = chem.Calibration(lambda time: 10 ** (-0.6 * time + 10.644), lb=160, ub=1_090_000, name="RI calibration")

    # load data and create signals
    file_name = r"C:\Users\nicep\Desktop\post_doc_2022\Data\Instrument\polymerization\RAFT9-L20-5min_N2_push.csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    df = df.drop(["Y:.10", "Y:.11", "Y:.12"], axis=1)
    signals = []
    for col in df.columns:
        signals.append(chem.SECSignal(name=col, ser=df[col], cal=cal_RI,
                                      x_label="retention time (min)", y_label="signal"))

    # process signals
    for sig in signals:
        sig.pipeline.add(chem_bc.adaptive_polynomial_baseline, remove_amount=0.7,)
        sig.peak_picking(lb=10, ub=13.1)

    # plotting
    colors = chem_pf.get_plot_color(len(signals))
    fig = signals[-1].plot(auto_open=False, op_cal=False, normalize=True, color=colors[0])
    op_cal = False
    for sig, color in zip(reversed(signals[:-1]), colors[1:]):
        if color == colors[-1]:
            op_cal = True
        fig = sig.plot(fig, auto_open=False, op_cal=op_cal, normalize=True, color=color)
    # fig.update_layout(legend=dict(x=.5, y=.95))
    # fig.data[4].fillcolor = 'rgba(172,24,25,0.5)'
    fig.update_xaxes(range=[8, 16])
    fig.write_html('temp.html', auto_open=True)

    # print out table
    df_table = signals[0].stats_df()
    for sig in signals[1:]:
        df_table = pd.concat([df_table, sig.stats_df()], axis=1)

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    print(df_table)


if __name__ == '__main__':
    main()
    print("done")
