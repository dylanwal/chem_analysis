
import pandas as pd
import well_plate

df = pd.read_csv(r"G:\Other computers\My Laptop\post_doc_2022\Data\Instrument\polymerization\RAFT_data\computed_data"
                 r".csv", index_col=0)

wp = well_plate.WellPlate(384, "rect")
wp.add_data(df["mw_n"])
wp.plot(key="mw_n")
