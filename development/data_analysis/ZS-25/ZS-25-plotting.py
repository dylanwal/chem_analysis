import numpy as np
import polars as pl
import plotly.graph_objects as go


def main():
    df = pl.read_csv(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\ZS-25\ZS-25_Diketone-slow-addition-t_mmols.csv")


    fig = go.Figure()
    for col in df.columns:
        if col == "TCB": continue
        if col == "time": continue
        fig.add_scatter(x=df["time"], y=df[col], name=col, mode="lines+markers")

    fig.show('browser')





if __name__ == "__main__":
    main()
