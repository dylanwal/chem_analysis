"""
testing how I calc stats.
"""
import plotly.graph_objs as go
import numpy as np
import scipy.stats as stats
from scipy.stats import kurtosis, skew, tstd


x = np.linspace(-15, 15, 50000)
distnames = ['laplace', 'norm', 'norm', 'norm', 'uniform', "skewnorm", "skewnorm"]

fig = go.Figure()

count = 0
count2=0
for distname in distnames:
    print(distname)
    if distname == 'uniform':
        dist = getattr(stats, distname)(loc=-2, scale=4)
    elif distname == "skewnorm":
        if count == 0:
            dist = getattr(stats, distname)(a=5)
            count += 1
        else:
            dist = getattr(stats, distname)(a=-3)
    elif distname == "norm":
        if count2 == 0:
            dist = getattr(stats, distname)
            count2 += 1
        elif count2==1:
            dist = getattr(stats, distname)(loc=1)
            count2 += 1
        else:
            dist = getattr(stats, distname)(scale=1.5)
    else:
        dist = getattr(stats, distname)

    data = dist.rvs(size=1_000_000)
    kur = kurtosis(data, fisher=True)
    s = skew(data)
    std1 = tstd(data)
    print(f"{kur}; {s}; {std1}")
    y = dist.pdf(x)

    mean = np.trapz(x=x, y=x * y)
    std = np.sqrt(np.trapz(x=x, y=y * (x - mean) ** 2))
    skews = np.trapz(x=x, y=y * (x - mean) ** 3) / std ** 3
    k = np.trapz(x=x, y=y * (x - mean) ** 4) / std ** 4 - 3
    print(f" {k}; {skews}; {std}; {mean}")


    fig.add_trace(go.Scatter(x=x, y=y, name="{}, {}".format(distname, round(kur, 3))))

fig.write_html("temp.html", auto_open=True)
