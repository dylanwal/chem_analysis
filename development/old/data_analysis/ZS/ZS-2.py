import chem_analysis as ca
import plotly.graph_objects as go

from development.plotly_utils import *


def main():
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(r"C:\Users\nicep\Desktop\research_wis\data\ZS\2\2\ZS-2-t1380.D")

    fig = go.Figure()
    fig.add_scatter(x=fid.x, y=fid.y)
    fig.layout.yaxis.range = (-10_000, 400_000)
    fig.layout.xaxis.range = (5, 40)
    fig.layout.xaxis.title.text = "<b>time (min)"
    fig.layout.yaxis.title.text = "<b>intensity"
    fig.show()


if __name__ == '__main__':
    main()
