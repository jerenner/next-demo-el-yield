"""Plot the quantities shown in Figures 2–5 of the preprint."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from .events import gaussian
from .yield_fit import relative_values

FIELD_LABEL = r"$E/p$ (kV/cm/bar)"


def save(figure, path):
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def yield_figure(fits, output):
    fig = plt.figure(figsize=(8, 6), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[4, 1], width_ratios=[1, 8])
    main = fig.add_subplot(grid[0, :])
    low = fig.add_subplot(grid[1, 0])
    high = fig.add_subplot(grid[1, 1])
    low.set(xlim=(0.87, 1.1), ylim=(-5, 25), ylabel=r"Residual / $\sigma_y$")
    high.set(xlim=(1.1, 3.08), ylim=(-2, 2), yticks=[-2, 0, 2])
    high.yaxis.tick_right()
    fig.supxlabel(FIELD_LABEL, fontsize=11)
    main.set(xlim=(0.87, 3.08), ylabel=r"$Y/(p d)$ (photons/e$^{-}$/cm/bar)")
    for fit, color in zip(fits, ("#46327e", "#277f8e", "#35b779"), strict=True):
        label = (f"{fit.pressure_bar:.1f} bar\n"
                 f"A = {fit.slope:.1f} ± {fit.slope_error:.1f}, "
                 f"X₀ = {fit.threshold:.3f} ± {fit.threshold_error:.3f}\n"
                 f"χ²/ndf = {fit.chi2_ndf:.2f}")
        main.errorbar(fit.field, fit.reduced_yield, xerr=fit.field_error, yerr=fit.yield_error,
                      fmt="o", ms=4, color=color, capsize=2, label=label)
        x = np.linspace(fit.field.min(), fit.field.max(), 100)
        main.plot(x, fit.predict(x), "--", color=color)
        residual = (fit.reduced_yield - fit.predict(fit.field)) / fit.yield_error
        for panel in (low, high):
            panel.plot(fit.field, residual, "o", ms=4, color=color)
    main.legend(loc="upper left")
    for panel in (low, high):
        panel.axhline(0, color="0.5", lw=0.8)
    for panel in (main, low, high):
        panel.grid(alpha=0.25)
    save(fig, output)


def relative_figure(fits, parameter, output):
    pressures = [f.pressure_bar for f in fits[1:]]
    values, errors = relative_values([getattr(f, parameter) for f in fits],
                                     [getattr(f, parameter + "_error") for f in fits])
    fig, ax = plt.subplots(figsize=(8, 3), layout="constrained")
    if parameter == "slope":
        # Freitas et al. (preprint reference 10), Table 1.
        literature, literature_err = relative_values(
            [141, 141, 142, 151, 161, 170, 162], [6, 7, 5, 5, 9, 10, 10]
        )
        ax.errorbar([3, 4, 5, 6, 8, 10], literature, yerr=literature_err,
                    fmt="ks", capsize=4, label="Freitas et al.")
        label = r"Relative slope $A/A_{2\,\mathrm{bar}}$"
    else:
        label = r"Relative threshold $X_0/X_{0,2\,\mathrm{bar}}$"
    ax.errorbar(pressures, values, yerr=errors, fmt="o", capsize=4, label="This work")
    ax.axhline(1, ls=":", color="0.4")
    ax.set(xlabel="Pressure (bar)", ylabel=label)
    ax.legend()
    ax.grid(alpha=0.25)
    save(fig, output)


def kr_example(candidates, fit, point, output_dir):
    near_peak = candidates[np.abs(candidates["S2e"] - fit.mean) < 3 * fit.width]
    fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
    _, _, _, counts = ax.hist2d(near_peak["X"], near_peak["Y"], bins=50,
                                range=((-100, 100), (-100, 100)), norm=LogNorm(vmin=1))
    fig.colorbar(counts, ax=ax, label="Counts/bin")
    ax.add_patch(plt.Circle((0, 0), 20, fill=False, ls="--", color="red", label="R < 20 mm"))
    ax.set(aspect="equal", xlabel="Reconstructed X (mm)", ylabel="Reconstructed Y (mm)")
    ax.legend()
    save(fig, output_dir / f"run_{point.run}_xy.pdf")

    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    ax.stairs(fit.counts, fit.edges, label="Spectrum")
    ax.errorbar(fit.fit_centers, fit.fit_counts, yerr=np.sqrt(fit.fit_counts),
                fmt="k.", ms=4, label="Fit bins")
    x = np.linspace(point.fit_min_pe, point.fit_max_pe, 200)
    # The original Fig. 3 overlay normalizes the Gaussian to the tallest fit bin.
    curve = gaussian(x, fit.area, fit.mean, fit.width)
    ax.plot(x, curve * max(fit.fit_counts) / max(curve), "r-", label="Gaussian fit")
    ax.text(0.98, 0.95, f"μ = {fit.mean:.1f} ± {fit.mean_error:.1f} p.e.\n"
                       f"χ²/ndf = {fit.chi2_ndf:.2f}", ha="right", va="top",
            transform=ax.transAxes)
    ax.set(xlabel="S2 energy (p.e.)", ylabel="Counts/bin")
    ax.legend(loc="upper left")
    save(fig, output_dir / f"run_{point.run}_s2_fit.pdf")
