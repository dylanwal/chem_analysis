import numpy as np
import polars as pl
import plotly.graph_objects as go


def main():
    df = pl.read_csv(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\DJW-11-115\DJW-11-115-V1-t_mmols.csv")

    fig = go.Figure()
    for col in df.columns:
        if col == "TCB": continue
        if col == "time": continue
        fig.add_scatter(x=df["time"], y=df[col], name=col, mode="lines+markers")

    fig.show('browser')

    # classes
    carboxylic_acids = [f"CA{i}" for i in [4, 5, 6, 7, 8, 10]]
    di_carboxylic_acids = [f"DA{i}" for i in range(4, 7)]
    df = df.with_columns([
        pl.sum_horizontal(carboxylic_acids).fill_null(0).alias("total_CA"),
        pl.sum_horizontal(di_carboxylic_acids).fill_null(0).alias("total_DA")
    ])

    keto_acids = ["KA5_2", "KA6_2", "KA6_3"]
    df = df.with_columns([
        pl.sum_horizontal(keto_acids).fill_null(0).alias("total_KA"),
    ])
    ketone = ["K10_2", "K10_3", "K10_4", "K10_5"]
    df = df.with_columns([
        pl.sum_horizontal(ketone).fill_null(0).alias("total_K"),
    ])

    fig = go.Figure()
    fig.add_scatter(x=df["time"], y=df["decane"], name="decane")
    fig.add_scatter(x=df["time"], y=df["total_CA"], name="total_CA")
    fig.add_scatter(x=df["time"], y=df["total_DA"], name="total_DA")
    fig.add_scatter(x=df["time"], y=df["total_KA"], name="total_KA")
    fig.add_scatter(x=df["time"], y=df["total_K"], name="total_K")
    fig.show('browser')

    df.write_csv(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\DJW-11-115\DJW-11-115-V1-final.csv")


if __name__ == "__main__":
    main()
