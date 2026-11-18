# Methods and provenance

## From waveforms to Table 2

The production chain used the NEXT Collaboration's Invisible Cities software. Irene restores the PMT baseline, applies run-selected PMT and SiPM calibrations, and finds the S1 and S2 peaks. Sophronia records the reconstructed event quantities in `/DST/Events`. The analysis here consumes those tables. The source analysis archive preserves the IC working copy at revision `315ed2e0fd8254db9af44f0815cf1e4da966327b`, its local changes, and the exact `localdb.DEMOPPDB.sqlite3`.

The calibration audit matched the embedded `/DB/DataPMT` and `/DB/DataSiPM` records in 2,692 Irene and 2,692 Sophronia files across the 83 paper runs. The three PMT gains (34.6135, 33.9873, and 41.4920 ADC/p.e.) come from calibration run 13727. The 256 SiPM gain entries (16.1258–18.8804 ADC/p.e.) come from run 13691. The database selects both gain sets from run 13727 through 100000, inclusive. One SiPM, sensor 14033, is masked. The calibration spectra, fit scripts, database, and file-by-file audit are in the legacy EL yield calibration archive.

The published peak reconstruction follows `scan_paper_v4/scan_all_pressures.ipynb` (zero-based cells 2, 3, 5, 6, and 9). It selects one S1 and one S2, keeps events with 100–32,000 p.e. in S2, and cuts on `0 < X² + Y² < 400 mm²` and `100 < Z < 300`. `Z` is the reconstructed drift coordinate used by the original selection. The separate drift-time studies use `DT`. The run-specific 80-bin histogram ranges and peak-fit windows are in `data/scan_points.csv`. The fit uses a Gaussian and weights each populated bin by its counting uncertainty, √count. No lifetime or geometric correction is applied to the S2 peak fits. The archived notebook's missing-map path labeled these uncorrected energies `S2ec_pes`.

The 83 saved means and errors came from the eight v4 `scan_Nbar.npz` products for nominal pressures 2–9 bar. Their rounded values and voltages match all 83 entries extracted directly from the arXiv v1 PDF's Table 2. Refitting the archived event tables reproduced the saved means and errors within 0.00000007 p.e. Figure 3 uses run 14017. The original plotting cells are 12 and 13 of `scan_all_pressures.ipynb`, and the manuscript's source PDFs match the archived v3 plots. The fitted mean is **7497.25 ± 7.12 p.e.** with χ²/ndf ≈ 1.25.

## Yield and pressure fits

The summary reproduces the v4 `scan_summary.ipynb` fit and ratio calculations (cells 2–5, 8, 12, 14, and 19). For a gate voltage in kV, a 1 cm gap, and a pressure in bar, the plotted field is `V/(p d)`. The reduced yield is `S2/(εLC Ne p d)` with `Ne = 41560/21.9`. The input pressure errors enter both axis uncertainties. The vertical error also combines the S2 statistical error and the 1.5% temporal systematic in quadrature.

Only points with `E/p > 1.1 kV/cm/bar` enter the fit. Figure 4 also shows the lower-field points in the residual panel. The code uses orthogonal distance regression for `A * (E/p) + B` and derives `X0 = -B/A`. The covariance is the unscaled ODR covariance, matching the original notebook. The plotted residuals use vertical errors. The figure uses the scans at 2.0052, 6.1571, and 9.4157 bar. The relative plots divide by the lowest-pressure fit, propagate its uncertainty, and omit that reference point from the plotted ratios. The comparison measurements are from Freitas et al., Table 1 (preprint reference 10).

The fit code calculates `εLC = 4.33 * 39.9 / 41560 ≈ 0.0041570`, which gives slopes of approximately 211.3 and 221.4 photons/e⁻/kV at the endpoint pressures, matching the figure labels. The arXiv v1 caption gives 0.0045. That value scales the absolute yields and slopes by about 0.924 and leaves the relative slopes and thresholds essentially unchanged. The global light-collection uncertainty is common to all points and is described separately in the preprint.

The saved weighted line fit of the seven relative slopes gives about 3.76 for slope divided by its fit error. The preprint reports 3.7σ. The implementation retains the notebook's use of pointwise ratio errors in this trend fit. The tabulated density interpolation at 293 K is from the v4 xenon density data (NIST Chemistry WebBook). Applying the ratio of pressure to density gives equivalent thresholds of about 0.916 and 0.848 at the endpoint pressures.

## Supporting calculations

The S1 study follows `scan_paper_v4/s1_extrapolation.ipynb`. For run 14046 it selects a krypton S2 energy band and uses the trend of squared S2 width with drift time to remove misidentified pulses. A narrower spatial sample and 100 iterations of a linear profile fit give an S1 intercept of about 4.33 p.e. at zero drift time. Its formal fit error is about 0.01 p.e. and the preprint uses 1.01 p.e. for the systematic uncertainty.

The electron lifetime diagnostic follows the preliminary linear band and final exponential profile fit in `scan_all_pressures.ipynb` cell 6. Individual fits can be poorly constrained when the S2 profile is nearly flat. The resulting values are diagnostics and are not used to correct the published yields. The separate COMSOL model and campaign slow-control records live with the original analysis material.

## Archive inputs

Use the following files relative to the original analysis directory. The archive retains the complete file hierarchy and calibrated production datasets.

| Step | Archive inputs |
| --- | --- |
| S2 refits and S1/lifetime studies | `<run>/sophronia/run_<run>_*_ldc1_trg0.v2.3.1.20250218.Kr.sophronia.h5`, `/DST/Events` |
| Figure 2, event 1886 | `13980/data/run_13980_0000_ldc1_trg0.waveforms.h5` and the corresponding file `13980/irene/run_13980_0000_ldc1_trg0.v2.3.1.20250218.Kr.irene.h5` |
| Figure 2, PMT calibration | The archived IC environment and its DEMOPP database, selected for run 13980 |

The waveform uses the historical IC baseline restoration (first 14,000 samples, threshold 5) and the database gains. The Irene peak tables supply the S1 and S2 time windows. The 25 ns samples for the displayed 320–480 µs segment are bundled in `data/figure_2_waveform.npz`. To recreate that file, run this command in the historical Linux IC environment (`IC-3.8-2022-04-13`, Python 3.8.19, with `ICTDIR` pointing at the archived IC working copy). Then pass the result to the modern `waveform` plotting command using `--waveform-data`.

```sh
python -m next_demo_el_yield.waveform --analysis-root "$ANALYSIS_ROOT" --output /tmp/figure_2_waveform.npz
```

The full detector inputs, original notebooks, calibration database, and operating logs are retained in the off-repository archive.
