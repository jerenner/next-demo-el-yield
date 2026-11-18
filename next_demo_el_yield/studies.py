"""Supporting studies of S1 light collection and electron attachment."""

import numpy as np
from scipy.optimize import curve_fit

from .events import fiducial, kr_candidates, load_events
from .yield_fit import line


def profile(x, y, limits, bins):
    edges = np.linspace(*limits, bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    groups = np.digitize(x, edges) - 1
    rows = []
    for i, center in enumerate(centers):
        values = y[groups == i]
        if len(values):
            error = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else np.nan
            rows.append((center, values.mean(), error))
    return np.array(rows).reshape(-1, 3).T


def weighted_line(profile_result):
    x, mean, error = profile_result
    valid = np.isfinite(error) & (error > 0)
    return curve_fit(line, x[valid], mean[valid], sigma=error[valid], absolute_sigma=True)


def extrapolate_s1(root):
    run = 14046
    events = load_events(root, run, ["nS1", "nS2", "S2e", "X", "Y", "Z", "DT", "S2w", "S1e"])
    events = kr_candidates(events)
    peak = (events["S2e"] >= 11800) & (events["S2e"] <= 14500)
    width_sample = events[peak & (events["DT"] < 280) &
                          (events["S2w"] >= 0) & (events["S2w"] <= np.sqrt(250))]
    t, width2 = width_sample["DT"], width_sample["S2w"]**2
    x, mean, _ = profile(t, width2, (0, 280), 25)
    preliminary = np.polyfit(x, mean, 1)
    residual = width2 - line(t, *preliminary)
    loose = np.abs(residual) < np.std(residual)
    width_line, _ = weighted_line(profile(t[loose], width2[loose], (0, 280), 25))
    width_spread = np.std(width2[loose] - line(t[loose], *width_line))

    events = fiducial(events[peak], radius_mm=40, z_interval=(0, 300))
    selected = (np.abs(events["S2w"]**2 - line(events["DT"], *width_line)) < 1.5 * width_spread)
    selected &= ((events["DT"] > 0) & (events["DT"] < 300)
                 & (events["S1e"] > 0) & (events["S1e"] < 15))
    t, s1 = events["DT"][selected], events["S1e"][selected]
    band = np.ones(t.shape, dtype=bool)
    for _ in range(100):
        params, covariance = weighted_line(profile(t[band], s1[band], (20, 300), 25))
        residual = s1 - line(t, *params)
        band = np.abs(residual) < np.std(residual)
    return dict(run=run, slope_pe_per_us=float(params[0]), intercept_pe=float(params[1]),
                intercept_fit_error_pe=float(np.sqrt(covariance[1, 1])),
                intercept_systematic_pe=1.01)


def electron_lifetime(root, point):
    events = kr_candidates(load_events(root, point.run, ["nS1", "nS2", "DT", "S2e"]))
    low, high = point.fit_min_pe, point.fit_max_pe
    t, energy = events["DT"], events["S2e"]
    time_window = (t >= 50) & (t <= 300)
    broad = time_window & (energy >= 0.75 * low) & (energy <= 1.25 * high)
    x, mean, _ = profile(t[broad], energy[broad], (50, 300), 10)
    preliminary = np.polyfit(x, mean, 1)
    band = time_window & (np.abs(energy - line(t, *preliminary)) < 0.6 * (high - low))
    x, mean, error = profile(t[band], energy[band], (50, 300), 10)
    seed = (mean[0], -x[-1] / np.log(mean[-1] / mean[0]))
    params, covariance = curve_fit(lambda t, e0, lifetime: e0 * np.exp(-t / lifetime),
                                   x, mean, p0=seed, sigma=error, maxfev=5000)
    return dict(run=point.run, lifetime_us=float(params[1]),
                lifetime_error_us=float(np.sqrt(covariance[1, 1])))
