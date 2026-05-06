import pathlib

import polars as pl
import polars.selectors as cs

rf = {
    "TCB": 6,
    "decane": 10,
    "CA3": 5,
    "CA4": 6,
    "CA5": 7,
    "CA6": 8,
    "CA7": 9,
    "CA8": 10,
    "CA9": 11,
    "CA10": 12,
    "DA4": 8,
    "DA4_double": 7.9,
    "DA5": 9,
    "DA6": 10,
    "DA7": 11,
    "DA8": 12,
    "DA9": 13,
    "DA10": 14,
    "K10_2": 9,
    "K10_3": 9,
    "K10_4": 9,
    "K10_5": 9,
    "A10_1": 12.5,
    "A10_2": 12.5,
    "A10_3": 12.5,
    "A10_4": 12.5,
    "A10_5": 12.5,
    "AA10_1": 10.5,
    "AA10_2": 10.5,
    "AA10_3": 10.5,
    "AA10_4": 10.5,
    "AA10_5": 10.5,
    "KA5_2": 6,
    "KA6_3": 7,
    "KA6_2": 7,
    "HA2": 6.5,
    "HAA2": 4.5,
    "HA3": 7.5,
    "HAA3": 5.5,
    "HA3_2": 7.5,
    "HAA3_2": 5.5,
    "L5_4": 4,
    "K9_2": 8,
    "K8_2": 7,
    "K7_2": 6,
    "K10_4,K10_5,KA5_2": 7,
    "DK10_2_1": 14.3,
    "DK10_2_2": 14.3,
    "DK10_2_3": 14.3,
    "DK10_1_1": 11.15,
    "DK10_1_2": 11.15,
}


def main():
    path = pathlib.Path(
        r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\ZS-25\ZS-25_Diketone-slow-addition-t_areas.csv")
    df_area = pl.read_csv(path)
    df_area = df_area.fill_null(0)
    TCB_conc = 0.053  # M
    volume = 0.000025
    multiplier = TCB_conc * volume * 15 / 0.000025
    print(df_area)

    print("Keys not in response factor:")
    for col in df_area.columns:
        if col not in rf.keys():
            print(col)
    df_mmol = df_area.with_columns(
        pl.col(col) / pl.col("TCB") * (rf["TCB"] / rf.get(col, 10)) * multiplier for col in df_area.columns if
        col != "time")

    df_mmol = (
        df_mmol
        .with_columns(
            # 1. Create the sum of the DK10 columns
            pl.sum_horizontal(cs.starts_with("DK10")).alias("tDK10_2_4")
        )
        .select(
            # 2. Keep everything EXCEPT the original DK10 columns
            ~cs.starts_with("DK10")
        )
    )

    print(df_mmol)
    df_mmol.write_csv(path.with_stem(path.stem.replace("area", "mmol")))


if __name__ == "__main__":
    main()
