def local_run():
    from scipy.stats import norm
    n = 1000
    rv = norm(loc=n/2, scale=10)
    x = np.linspace(0, n, n)
    y = np.linspace(0, n, n) + 20 * np.random.random(n) + 5000 * rv.pdf(x)
    signal = Signal(name="test", x=x, y=y, x_label="x_test", y_label="y_test")
    signal.baseline(deg=1)
    signal.auto_peak_picking()
    signal.stats()
    signal.plot()
    print("done")


if __name__ == '__main__':
    local_run()