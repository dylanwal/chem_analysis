import pathlib


import chem_analysis as ca
import plotly.graph_objs as go


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / "GCMS" / f"{label}_data.ms"
    fid_file = data_path / "GCMS" / f"{label}_FID1A.ch"
    ini_file = data_path / "GCMS" / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


ms, fid = load_single_file(pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_81"), "DJW-11-81-V2-t60-PPh3")

fig = go.Figure()
fig.add_scatter(x=fid.x, y=fid.y, mode="lines")
fig.layout.showlegend =False
fig.layout.xaxis.title.text = "<b>time (min)</b>"
fig.layout.yaxis.title.text = "<b>signal</b>"
fig.layout.xaxis.title.font = dict(size=24)
fig.layout.xaxis.tickfont = dict(size=18)
fig.layout.yaxis.title.font = dict(size=24)
fig.layout.yaxis.tickfont = dict(size=18)
fig.layout.xaxis.range = [4, 47]
fig.layout.yaxis.range = [-100_000, 450_000]
fig.layout.height = 600
fig.layout.width = 1000

fig.show()