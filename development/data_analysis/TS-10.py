import pathlib

import numpy as np
import polars as pl

import chem_analysis as ca


folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\TS\TS-11\GCMS")
data_label = "TS-11-t" #"TS7-d4-ketone-t"
times = np.array([10, 20, 30, 40, 50, 60, 75, 90, 105, 120, 180, 240])  #


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

def main():
    labels = [data_label + str(i) for i in times]
    raw_data: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in labels ]
    print(raw_data)

    ms_spectras = []
    for i in range(len(raw_data)):
        ms = raw_data[i][0].extract_ms(30.2,30.7)
        ms = ms.to_numpy()
        index = np.argmin(np.abs(ms[:,0] - 180))
        ms_spectras.append(ms[:index])

    # print(all(i.shape[0] == ms_spectras[0].shape[0] for i in ms_spectras))

    df = pl.DataFrame({"mz":ms_spectras[0][:,0]} | {str(i): ms_spectras[i][:,1] for i in range(len(ms_spectras))})

    import polars.selectors as cs

    filtered_df = df.filter(
        ~pl.all_horizontal(cs.all() - cs.first() == 0)
    )
    filtered_df.write_csv("TS11-ms.csv")
    print(filtered_df)
    # fig = ca.plot.signal(ms)
    # fig.show("browser")


if __name__ == "__main__":
    main()