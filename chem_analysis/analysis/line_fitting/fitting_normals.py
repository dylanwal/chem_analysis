from typing import Callable, Sequence, Any
import logging

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm, qmc
from sklearn.metrics import r2_score

logger = logging.getLogger(__name__)


# Input || x: np.ndarray, args: Sequence
# Return || y: np.ndarray
FittingFunction = Callable[[np.ndarray, Sequence[Any]], np.ndarray]


def n_normals(x, *params):
    n = len(params) // 3  # Each normal has 3 params: amplitude, mean, std
    y = np.zeros_like(x)
    for i in range(n):
        amplitude = params[i * 3]
        mu = params[i * 3 + 1]
        sigma = params[i * 3 + 2]
        y += amplitude * norm.pdf(x, loc=mu, scale=sigma)
    return y


def fit_n_normals(num_normals: int, x_data: np.ndarray, y_data: np.ndarray, num_trials=5) -> tuple | None:
    """Fit models with 1, 2, or 3 normal distributions"""
    best_fit = None
    best_bic = np.inf

    # Latin Hypercube Sampling for mean positions
    sampler = qmc.LatinHypercube(d=num_normals)
    sample = sampler.random(num_trials)
    mean_bounds = (x_data.min(), x_data.max())
    sampled_means = qmc.scale(sample, mean_bounds[0], mean_bounds[1])

    for trial in range(num_trials):
        # Generate initial guesses using sampled means
        p0 = []
        for i in range(num_normals):
            p0 += [1/num_normals, sampled_means[trial, i], x_data.std()]  # amplitude, mean, std

        try:
            popt, _ = curve_fit(n_normals, x_data, y_data, p0=p0)
            # Calculate residuals and BIC
            residuals = y_data - n_normals(x_data, *popt)
            sse = np.sum(residuals**2)
            k = len(popt)  # Number of parameters
            n = len(x_data)
            bic = n * np.log(sse / n) + k * np.log(n)  # BIC formula

            if bic < best_bic:
                best_bic = bic
                best_fit = (popt, bic)

        except RuntimeError:
            continue  # Skip failed fits

    if best_fit is None:
        return None

    r2 = r2_score(y_data, n_normals(x_data, *best_fit[0]))
    if r2 < 0.9:
        logger.info('fit was not good enough')
        return None

    return best_fit


def integrate_by_fitting_single(
    x: np.ndarray,
    y: np.ndarray,
    num_normals: tuple[int, int] = (1, 3),
):
    # fit models
    models = {}
    for num_normals_ in range(*num_normals):
        fit_result = fit_n_normals(num_normals_, x, y, num_trials=20)
        if fit_result:
            models[num_normals] = fit_result

    # Select the best model based on Bayesian Information Criterion (BIC)
    bic_values = {num_normals: fit[1] for num_normals, fit in models.items()}
    best_num_normals = min(bic_values, key=bic_values.get)
    best_popt = models[best_num_normals][0]
    logger.debug(f"Best fit: {best_num_normals} normal distributions (BIC: {bic_values[best_num_normals]:.2f})")

    r2 = r2_score(y, n_normals(x, *best_popt))
    if r2 < 0.9:
        logger.info('fit was not good enough')
        return None

    # Generate the fit for the best model
    #y_fit = n_normals(x, *best_popt)

    return best_popt
