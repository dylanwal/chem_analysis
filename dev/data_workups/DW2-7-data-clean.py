import pathlib

import numpy as np

import chem_analysis as ca

remove = [
    [-338, 0],
    [5821, 6183],
    [12_824, 14_308],
    [20_960, 21_780],
    [29_272, 29_325],
    [30_584, 30_728],
]

path = pathlib.Path(r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\sec.csv")
data = np.loadtxt(path, delimiter=",")
print(data.shape)

for r in reversed(remove):  # start from back to front since we are shifting times
    slice_ = ca.utils.general_math.get_slice(data[:, 0], *r)
    data = np.delete(data, slice_, axis=0)
    data[slice_.start:, 0] -= (r[1] - r[0])  # shift time so gap doesn't exist

print(data.shape)
path = path.with_stem(path.stem + "_proc")
np.savetxt(path, data, delimiter=",")
