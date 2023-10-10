
import numpy as np

import chem_analysis as ca


def cal_func(time: np.ndarray):
    return 10 ** (0.0167 * time ** 2 - 0.9225 * time + 14.087)


calibration = ca.sec.ConventionalCalibration(cal_func, lower_bound_mw=900, upper_bound_mw=319_000)


# load test data (log-normal distribution)
D = 1.05
Mn = 15_000  # g/mol
x = np.linespace(0, 30, 1_000)
y = 1 / (x * np.sqrt(2 * np.pi * np.log(D))) * np.exp(-1 * (np.log(x / Mn) + np.log(D) / 2) ** 2 / (2 * np.log(D)))

signal = ca.sec.SECSignal(x, y, calibration)

signal.print_stats()
