"""Reproduce the Kr S2 peak fit from reconstructed Sophronia event tables."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import curve_fit


def load_events(root, run, columns):
    pattern = f"run_{run}_*_ldc1_trg0.v2.3.1.20250218.Kr.sophronia.h5"
    files = sorted((Path(root) / str(run) / "sophronia").glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} under {root}/{run}/sophronia")
    arrays = []
    for path in files:
        with h5py.File(path, "r") as handle:
            arrays.append(handle["DST/Events"].fields(columns)[:])
    return np.concatenate(arrays)


def kr_candidates(events):
    return events[(events["nS1"] == 1) & (events["nS2"] == 1)
                  & (events["S2e"] >= 100) & (events["S2e"] <= 32000)]


def fiducial(events, radius_mm=20, z_interval=(100, 300)):
    radius_squared = events["X"]**2 + events["Y"]**2
    return events[(radius_squared > 0) & (radius_squared < radius_mm**2)
                  & (events["Z"] > z_interval[0]) & (events["Z"] < z_interval[1])]


def gaussian(x, area, mean, sigma):
    if sigma <= 0:
        return np.inf
    return area / (2 * np.pi)**0.5 / sigma * np.exp(-0.5 * (x - mean)**2 / sigma**2)


@dataclass(frozen=True)
class PeakFit:
    counts: np.ndarray
    edges: np.ndarray
    fit_centers: np.ndarray
    fit_counts: np.ndarray
    area: float
    mean: float
    width: float
    mean_error: float
    chi2_ndf: float


def fit_peak(energies, point):
    counts, edges = np.histogram(energies, bins=80, range=(point.hist_min_pe, point.hist_max_pe))
    centers = (edges[:-1] + edges[1:]) / 2
    low, high = point.fit_min_pe, point.fit_max_pe
    in_window = (centers >= low) & (centers <= high) & (counts > 0)
    x, y = centers[in_window], counts[in_window]
    if len(x) <= 3:
        raise ValueError(f"Too few populated bins to fit run {point.run}")
    inside = energies[(energies >= low) & (energies <= high)]
    width = 0.5 * np.std(inside) if len(inside) > 1 else (high - low) / 6
    if not np.isfinite(width) or width <= 0:
        width = (high - low) / 10
    seed = (max(y) * width * np.sqrt(2 * np.pi), x[np.argmax(y)], width)
    values, covariance = curve_fit(gaussian, x, y, p0=seed, sigma=np.sqrt(y), maxfev=10000)
    chi2 = np.sum(((y - gaussian(x, *values)) / np.sqrt(y))**2)
    return PeakFit(counts, edges, x, y, *values, np.sqrt(covariance[1, 1]), chi2 / (len(x) - 3))


def analyze_run(root, point):
    events = load_events(root, point.run, ["nS1", "nS2", "X", "Y", "Z", "S2e"])
    candidates = kr_candidates(events)
    selected = fiducial(candidates)
    # The paper peaks use uncorrected S2. A missing correction map in the
    # historical notebook fell back to this column under the name S2ec_pes.
    return candidates, fit_peak(selected["S2e"], point)
