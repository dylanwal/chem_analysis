import pathlib

import polars as pl

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
    "DA4_double":7.9,
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
}


ratio = [.21, .11, .67]
rf["K10_4,K10_5,KA5_2"] = rf.get("K10_4") * ratio[0] + rf.get("K10_5") * ratio[1] + rf.get("KA5_2") * ratio[2]


def main():
    path = pathlib.Path(r"C:\Users\nicep\Desktop\pyth_proj\chem_analysis\development\data_analysis\DJW-11-115\DJW-11-115-V1-t_areas.csv")
    df_area = pl.read_csv(path)
    df_area = df_area.fill_null(0)
    TCB_conc = 0.053 #M
    volume = 0.000025
    multiplier = TCB_conc * volume * 15/0.000025
    print(df_area)


    print("Keys not in response factor:")
    for col in df_area.columns:
        if col not in rf.keys():
            print(col)
    df_mmol = df_area.with_columns(pl.col(col)/pl.col("TCB")*(rf["TCB"]/rf.get(col, 10))*multiplier for col in df_area.columns if col !="time")


    # Create K5 (90%) and k6 (10%)
    df_mmol = df_mmol.with_columns([
        (pl.col("K10_4,K10_5,KA5_2") * ratio[0]).alias("K10_4"),
        (pl.col("K10_4,K10_5,KA5_2") * ratio[1]).alias("K10_5"),
        (pl.col("K10_4,K10_5,KA5_2") * ratio[2]).alias("KA5_2"),
    ])

    # Optional: If you want to remove the original 'kts' column
    df_mmol = df_mmol.drop("K10_4,K10_5,KA5_2")


    print(df_mmol)
    df_mmol.write_csv(path.with_stem(path.stem.replace("area", "mmol")))


if __name__ == "__main__":
    main()
