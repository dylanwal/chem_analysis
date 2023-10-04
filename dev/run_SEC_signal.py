
def local_run():
    from scipy.stats import norm

    def cal_func(time: np.ndarray):
        return 10**(0.0167 * time ** 2 - 0.9225 * time + 14.087)

    cal = Cal(cal_func, lb=900, ub=319_000)

    nx = 1000
    x = np.linspace(0, 25, nx)
    rv = norm(loc=15, scale=0.6)
    rv2 = norm(loc=18, scale=0.6)
    y = 5 * np.linspace(0, 1, nx) + np.random.random(nx) + 100 * rv.pdf(x) + 20 * rv2.pdf(x)

    sig = SECSignal(name="RI data", x=x, y=y, x_label="time (min)", y_label="intensity", cal=cal)
    sig.auto_peak_baseline(deg=1)
    sig.plot()
    sig.stats()
    print("done")



def local_run():
    from scipy.stats import norm
    import numpy as np

    def cal_func(time: np.ndarray):
        return 10**(0.0167 * time ** 2 - 0.9225 * time + 14.087)

    cal = Cal(cal_func, lb=900, ub=319_000)

    nx = 1000
    ny = 3
    x = np.linspace(0, 25, nx)
    y = np.empty((ny, nx))
    for i in range(ny):
        rv = norm(loc=15, scale=0.6)
        rv2 = norm(loc=18, scale=0.6)
        y[i, :] = 5 * np.linspace(0, 1, nx) + np.random.random(nx) + 100 * rv.pdf(x) + 20 * rv2.pdf(x)
    df = pd.DataFrame(data=y.T, index=x)
    df.columns = ["RI", "UV", "LS"]
    df.index.names = ["time"]

    chro = SECChrom(data=df, cal=cal)
    chro.auto_peak_baseline(deg=1)
    chro.plot()
    chro.stats()
    print("done")


if __name__ == "__main__":
    local_run()