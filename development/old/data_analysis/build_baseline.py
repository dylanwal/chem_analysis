import pathlib

import chem_analysis as ca
import numpy as np


def build_baseline(data_path: pathlib.Path, label: str):
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    baseline = ca.p.baseline.AdaptiveAsymmetricLeastSquared(save_result=True)
    fid.processor.add(
        ca.p.edit.CutSpans(((1.2, 11), (29, 31)), invert=True),
        ca.p.resampling.EveryN(step=[int(fid.x_raw.size/6000)]),
        ca.p.smoothing.Gaussian(sigma=10),
        baseline
    )

    fig = ca.plot.signal(fid, raw=True)
    fig = ca.plot.baseline(fid, fig=fig)
    min_, max_ = np.min(baseline.baseline), np.max(baseline.baseline)
    fig.layout.yaxis.range = min_ - (max_ - min_) / 4, max_ + (max_ - min_) / 2

    return baseline.x, baseline.baseline, fig


if __name__ == '__main__':
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_53\GCMS")
    label = "DJW-11-53-run1-t0"
    x, y, fig = build_baseline(data_path, label)
    fig.show()
