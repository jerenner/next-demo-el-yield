"""Commands for reproducing the NEXT-DEMO++ EL yield analysis."""

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from .yield_fit import DATA, FIGURE_EFFICIENCY, density_thresholds, fit_all
from .yield_fit import read_scan_points, relative_slope_trend


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(args):
    from .figures import relative_figure, yield_figure

    points = read_scan_points(args.scan_points)
    fits = fit_all(points, args.efficiency)
    expected_pressures = [2.0052, 3.0690, 4.0622, 5.1604, 6.1571, 7.2177, 8.2689, 9.4157]
    if [fit.pressure_bar for fit in fits] != expected_pressures:
        raise ValueError("The paper summary requires all eight published pressures")
    fields = ("pressure_bar", "slope", "slope_error", "threshold", "threshold_error", "chi2_ndf")
    fit_rows = [{field: getattr(fit, field) for field in fields} for fit in fits]
    if args.check:
        expected = np.genfromtxt(DATA / "reference_fits.csv", delimiter=",", names=True)
        for field in fields:
            np.testing.assert_allclose([row[field] for row in fit_rows], expected[field],
                                       rtol=5e-5, atol=2e-6, err_msg=field)
        print("All eight pressure fits match the archived figure calculations.")
    for row, fit, density_threshold in zip(fit_rows, fits, density_thresholds(fits), strict=True):
        row.update(fit_points=int(fit.included.sum()), density_scaled_threshold=density_threshold)
    write_csv(args.output / "yield_fits.csv", fit_rows)
    table = [{"pressure_bar": p.pressure_bar, "gate_voltage_v": round(p.gate_voltage_v),
              "peak_mean_pe": round(p.peak_mean_pe), "peak_error_pe": round(p.peak_error_pe)}
             for p in points]
    write_csv(args.output / "table_2.csv", table)
    trend = relative_slope_trend(fits)
    trend["light_collection_efficiency"] = args.efficiency
    (args.output / "relative_slope_trend.json").write_text(json.dumps(trend, indent=2) + "\n")
    yield_figure([fits[0], fits[4], fits[7]], args.output / "figure_4_yield.pdf")
    for parameter in ("slope", "threshold"):
        relative_figure(fits, parameter, args.output / f"figure_5_relative_{parameter}.pdf")


def refit(args):
    from .events import analyze_run
    from .figures import kr_example

    points = read_scan_points(args.scan_points)
    selected = set(args.run) if args.run else {p.run for p in points}
    missing = selected - {p.run for p in points}
    if missing:
        raise ValueError(f"Unknown runs: {sorted(missing)}")
    output = []
    for point in points:
        if point.run not in selected:
            continue
        candidates, fit = analyze_run(args.analysis_root, point)
        if args.check:
            np.testing.assert_allclose([fit.mean, fit.mean_error],
                                       [point.peak_mean_pe, point.peak_error_pe],
                                       rtol=0, atol=1e-3, err_msg=f"run {point.run}")
        output.append(asdict(replace(point, peak_mean_pe=fit.mean, peak_error_pe=fit.mean_error)))
        if point.run == 14017 or len(selected) == 1:
            kr_example(candidates, fit, point, args.output)
        print(f"Run {point.run}  mean {fit.mean:.3f} ± {fit.mean_error:.3f} p.e.", flush=True)
    write_csv(args.output / "refit_scan_points.csv", output)


def studies(args):
    from .studies import electron_lifetime, extrapolate_s1

    s1 = extrapolate_s1(args.analysis_root)
    (args.output / "s1_extrapolation.json").write_text(json.dumps(s1, indent=2) + "\n")
    points = read_scan_points(DATA / "scan_points.csv")
    write_csv(args.output / "electron_lifetimes.csv",
              [electron_lifetime(args.analysis_root, point) for point in points])


def waveform(args):
    from .waveform import plot_waveform

    plot_waveform(args.waveform_data, args.output / "figure_2_waveform.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True, dest="command")
    summary = commands.add_parser("summary", help="Reproduce Table 2 and Figures 4–5")
    peaks = commands.add_parser("refit", help="Refit reconstructed S2 peaks and reproduce Figure 3")
    diagnostic = commands.add_parser("studies", help="Run the S1 intercept and lifetime studies")
    example = commands.add_parser("waveform", help="Plot the calibrated Figure 2 waveform")
    for command in (summary, peaks, diagnostic, example):
        command.add_argument("--output", type=Path, default=Path("results"))
    for command in (peaks, diagnostic):
        command.add_argument("--analysis-root", type=Path, required=True)
    example.add_argument("--waveform-data", type=Path, default=DATA / "figure_2_waveform.npz")
    for command in (summary, peaks):
        command.add_argument("--scan-points", type=Path, default=DATA / "scan_points.csv")
        command.add_argument("--check", action="store_true", help="Compare against the saved analysis")
    summary.add_argument("--efficiency", type=float, default=FIGURE_EFFICIENCY)
    peaks.add_argument("--run", type=int, action="append", help="Refit this run, repeat for more runs")
    args = parser.parse_args()
    if args.command == "summary" and (not np.isfinite(args.efficiency) or args.efficiency <= 0):
        parser.error("--efficiency must be finite and positive")
    args.output.mkdir(parents=True, exist_ok=True)
    {"summary": summarize, "refit": refit, "studies": studies, "waveform": waveform}[args.command](args)
    print(f"Outputs written to {args.output}")


if __name__ == "__main__":
    main()
