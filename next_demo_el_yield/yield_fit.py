"""Convert Kr peak means to reduced yield and fit each pressure scan."""

import csv
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy.odr import Model, ODR, RealData
from scipy.optimize import curve_fit

DATA = Path(__file__).resolve().parents[1] / "data"
KR_ENERGY_EV = 41560.0
IONIZATION_WORK_EV = 21.9
S1_AT_ZERO_PE = 4.33
SCINTILLATION_WORK_EV = 39.9
FIGURE_EFFICIENCY = S1_AT_ZERO_PE * SCINTILLATION_WORK_EV / KR_ENERGY_EV
GAP_CM = 1.0
YIELD_SYSTEMATIC = 0.015
FIELD_CUT = 1.1


@dataclass(frozen=True)
class ScanPoint:
    run: int
    pressure_bar: float
    pressure_error_bar: float
    gate_voltage_v: float
    peak_mean_pe: float
    peak_error_pe: float
    hist_min_pe: float
    hist_max_pe: float
    fit_min_pe: float
    fit_max_pe: float


@dataclass(frozen=True)
class YieldFit:
    pressure_bar: float
    field: np.ndarray
    field_error: np.ndarray
    reduced_yield: np.ndarray
    yield_error: np.ndarray
    included: np.ndarray
    slope: float
    intercept: float
    slope_error: float
    threshold: float
    threshold_error: float
    chi2_ndf: float

    def predict(self, field):
        return self.slope * field + self.intercept


def read_scan_points(path):
    with Path(path).open(newline="") as stream:
        return [
            ScanPoint(**{f.name: f.type(row[f.name]) for f in fields(ScanPoint)})
            for row in csv.DictReader(stream)
        ]


def line(x, slope, intercept):
    return slope * x + intercept


def fit_scan(points, efficiency=FIGURE_EFFICIENCY):
    if not np.isfinite(efficiency) or efficiency <= 0:
        raise ValueError("Light-collection efficiency must be finite and positive")
    pressure = points[0].pressure_bar
    relative_pressure_error = points[0].pressure_error_bar / pressure
    field = np.array([p.gate_voltage_v for p in points]) / (1000 * GAP_CM * pressure)
    peak = np.array([p.peak_mean_pe for p in points])
    peak_error = np.array([p.peak_error_pe for p in points])

    # Follow the order of operations of the figure notebook.
    light_yield = peak / efficiency / (KR_ENERGY_EV / IONIZATION_WORK_EV)
    reduced_yield = light_yield / pressure / GAP_CM
    statistical = light_yield * (peak_error / peak)
    temporal_systematic = light_yield * YIELD_SYSTEMATIC
    combined = (statistical**2 + temporal_systematic**2)**0.5 / pressure / GAP_CM
    yield_error = (combined**2 + (reduced_yield * relative_pressure_error)**2)**0.5
    field_error = field * relative_pressure_error

    included = field > FIELD_CUT
    if included.sum() < 3:
        raise ValueError(f"Insufficient data to fit the scan at {pressure} bar")
    data = RealData(field[included], reduced_yield[included],
                    sx=field_error[included], sy=yield_error[included])
    fit = ODR(data, Model(lambda beta, x: line(x, *beta)), beta0=(140, -140 * 0.83))
    fit.set_job(fit_type=0)
    result = fit.run()
    if result.info not in (1, 2, 3):
        raise RuntimeError(f"ODR failed at {pressure} bar: {result.stopreason}")
    slope, intercept = result.beta
    covariance = result.cov_beta  # Historical notebook uses the unscaled ODR covariance.
    threshold = -intercept / slope
    jacobian = np.array([intercept / slope**2, -1 / slope])
    return YieldFit(pressure, field, field_error, reduced_yield, yield_error, included,
                    slope, intercept, np.sqrt(covariance[0, 0]), threshold,
                    np.sqrt(jacobian @ covariance @ jacobian), result.res_var)


def fit_all(points, efficiency=FIGURE_EFFICIENCY):
    grouped = {}
    for point in sorted(points, key=lambda p: (p.pressure_bar, p.gate_voltage_v)):
        grouped.setdefault(point.pressure_bar, []).append(point)
    return [fit_scan(group, efficiency) for group in grouped.values()]


def relative_values(values, errors):
    values, errors = np.asarray(values), np.asarray(errors)
    ratio = values / values[0]
    error = ratio * np.sqrt((errors / values)**2 + (errors[0] / values[0])**2)
    return ratio[1:], error[1:]


def relative_slope_trend(fits):
    """The saved notebook's weighted line, using the seven non-reference scans."""
    values, errors = relative_values([f.slope for f in fits], [f.slope_error for f in fits])
    pressure = [f.pressure_bar for f in fits[1:]]
    params, covariance = curve_fit(line, pressure, values, sigma=errors, absolute_sigma=True)
    slope_error = np.sqrt(covariance[0, 0])
    return dict(slope_per_bar=float(params[0]), slope_error_per_bar=float(slope_error),
                slope_over_error=float(params[0] / slope_error))


def density_thresholds(fits):
    """Express X0 in equivalent bar at the 293 K density of 1 bar."""
    density = np.genfromtxt(DATA / "xenon_density_293k.csv", delimiter=",", names=True)
    pressure, rho = density["pressure_bar"], density["density_kg_m3"]
    return [float(f.threshold * f.pressure_bar * np.interp(1, pressure, rho)
                  / np.interp(f.pressure_bar, pressure, rho)) for f in fits]
