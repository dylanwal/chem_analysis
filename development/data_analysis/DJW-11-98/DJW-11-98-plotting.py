import numpy as np
import polars as pl
import plotly.graph_objects as go


def main():
    df = pl.read_csv(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\DJW-11-98\DJW-11-98-t_mmols.csv")

    df = df.with_columns(
        pl.sum_horizontal(df.columns[1:]).fill_null(0).alias("mass balance"),
    )

    print(df)
    fig = go.Figure()
    for col in df.columns:
        if col == "TCB": continue
        if col == "time": continue
        fig.add_scatter(x=df["time"], y=df[col], name=col, mode="lines+markers")

    fig.show('browser')



if __name__ == "__main__":
    main()
