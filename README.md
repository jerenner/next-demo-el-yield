# NEXT-DEMO++ electroluminescence yield analysis

This repository documents and reproduces the analysis in J. Renner et al. (NEXT Collaboration), *Electroluminescence Yield Measurements in Xenon Gas with the NEXT-DEMO++ Detector*, [arXiv:2604.15905v1](https://arxiv.org/abs/2604.15905v1).

The analysis starts from run-dependent krypton peak fits in reconstructed Invisible Cities data and ends with the electroluminescence yield and its pressure dependence. The repository includes the 83 saved peak measurements, the run and fit-window settings, and the plotted samples for the example waveform. Use these inputs to regenerate Table 2 and Figures 2, 4, and 5 on a laptop. Figure 3 and the per-run peak refits use the archived detector data.

## Run from a clone

Use Python 3.11 or newer with the pinned dependencies. Run the commands from the repository root.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m next_demo_el_yield summary --check
python -m next_demo_el_yield waveform
```

The commands write to `results/`. The summary writes the rounded values of Table 2, all eight pressure fits with uncertainties, the relative-slope trend, and PDFs of the yield and relative fit parameters. The plot code redraws the published quantities in a compact layout. The `--check` option compares the numerical fits with the archived figure calculation.

The default light-collection efficiency is **0.0041570**, the value that reproduces the published figure slopes. The preprint text and Table 2 caption give **0.0045**. To follow the caption's normalization, run `python -m next_demo_el_yield summary --efficiency 0.0045 --output results/caption-normalization`. Both settings use the same peak measurements. See [the methods and provenance notes](docs/methods.md) for the effect on the absolute slope and for the provenance of the plots.

## Refit the detector data

Set `ANALYSIS_ROOT` to the directory containing numeric run directories, such as `13980/` and `14017/`. In the backed-up archive this is `pokopandegi/analysis/`.

```sh
python -m next_demo_el_yield refit --analysis-root "$ANALYSIS_ROOT" --run 14017 --check
python -m next_demo_el_yield refit --analysis-root "$ANALYSIS_ROOT" --check --output results/all-runs
python -m next_demo_el_yield summary --scan-points results/all-runs/refit_scan_points.csv --check
python -m next_demo_el_yield studies --analysis-root "$ANALYSIS_ROOT"
```

The first command draws the Figure 3 XY distribution and S2 fit. The second refits all 83 inputs from the Sophronia event tables. The studies command repeats the S1 light-collection extrapolation on run 14046 and calculates the S2-versus-drift-time lifetime diagnostics. See [the archive layout and waveform export instructions](docs/methods.md#archive-inputs).

## Where to start in the code

| File | Analysis step |
| --- | --- |
| [data/scan_points.csv](data/scan_points.csv) | Run list, gate voltages, pressures and errors, S2 fit windows, and saved peak results |
| [events.py](next_demo_el_yield/events.py) | Krypton selection, fiducial cuts, and Gaussian S2 fits |
| [yield_fit.py](next_demo_el_yield/yield_fit.py) | Yield conversion, orthogonal distance regression, and relative pressure trends |
| [studies.py](next_demo_el_yield/studies.py) | Supporting S1 extrapolation and electron lifetime fits |
| [waveform.py](next_demo_el_yield/waveform.py) | Figure 2 baseline restoration, PMT calibration, and peak windows |
| [figures.py](next_demo_el_yield/figures.py) | Figures 3–5 |
| [docs/methods.md](docs/methods.md) | Source notebooks, calibration provenance, inputs, and numerical cross-checks |
