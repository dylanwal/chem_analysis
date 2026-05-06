import numpy as np
import polars as pl
import plotly.graph_objects as go


def main():
    df = pl.read_csv(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\ZS-7-V2\ZS-7-V2-t_mmols.csv")


    fig = go.Figure()
    for col in df.columns:
        if col == "TCB": continue
        if col == "time": continue
        fig.add_scatter(x=df["time"], y=df[col], name=col, mode="lines+markers")

    fig.show('browser')


    fig = go.Figure()

    for col in df.columns:
        if col in ["TCB", "time"]:
            continue

        # 1. Get the initial concentration [M]0 (first data point)
        m0 = df[col][0]

        # 2. Safety Check: We can't divide by 0.
        # This skips products that start at 0 concentration.
        if m0 <= 0:
            continue

        # 3. Calculate ln([M]t / [M]0)
        # We use polars expressions for speed, or simple series math
        # formula: ln( current_val / initial_val )
        y_values = (df[col] / m0).log()

        fig.add_scatter(
            x=df["time"],
            y=y_values,
            name=col,
            mode="lines+markers"
        )

    # Add axis titles for clarity
    fig.update_layout(
        title="First Order Kinetic Plot",
        xaxis_title="Time",
        yaxis_title="ln([M]t / [M]0)"
    )

    fig.show('browser')


if __name__ == "__main__":
    main()
