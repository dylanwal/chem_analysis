
import numpy as np
import pandas as pd


def cal_Mn_D_from_wi(mw_i: np.ndarray, wi: np.ndarray) -> tuple[float, float]:
    """ calculate Mn and D from wi vs MW data (MW goes low to high) """
    data_points = len(mw_i)

    # flip data if giving backwards
    if mw_i[1] > mw_i[-1]:
        mw_i = np.flip(mw_i)
        wi = np.flip(wi)

    wi_d_mi = np.zeros(data_points)
    wi_m_mi = np.zeros(data_points)
    for i in range(data_points):
        if mw_i[i] != 0:
            wi_d_mi[i] = wi[i] / mw_i[i]
        wi_m_mi[i] = wi[i] * mw_i[i]

    mw_n = np.sum(wi) / np.sum(wi_d_mi)
    mw_w = np.sum(wi_m_mi) / np.sum(wi)
    mw_d = mw_w / mw_n
    return mw_n, mw_d


def main():
    file = r"C:\Users\nicep\Desktop\data.csv"
    df = pd.read_csv(file, header=None, index_col=0)

    range_ = (700, 80_000)
    lower = min(len(df.index), np.where(df.index < range_[0])[0][0])
    upper = max(0, np.where(df.index > range_[1])[0][-1])
    df_data = df.iloc[upper:lower]

    mw_i = df_data.index.to_numpy()

    wi = df_data[df_data.columns[0]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

    wi = df_data[df_data.columns[1]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

    wi = df_data[df_data.columns[2]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

    wi = df_data[df_data.columns[3]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

    wi = df_data[df_data.columns[4]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

    wi = df_data[df_data.columns[5]].to_numpy()
    wi2 = np.flip(wi / np.trapz(x=mw_i, y=wi))
    print(cal_Mn_D_from_wi(np.flip(mw_i), wi2))

if __name__ == "__main__":
    main()
