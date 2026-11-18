"""Figure 2 from one archived raw waveform and its Irene peak windows."""

import argparse
from pathlib import Path

import numpy as np

RUN = 13980
EVENT = 1886
SAMPLE_US = 0.025


def calibrate_waveform(root):
    """Use the historical IC baseline restoration and run-specific PMT gains."""
    import tables
    from invisible_cities.database.load_db import DataPMT
    from invisible_cities.evm.ic_containers import CalibVectors, DeconvParams
    from invisible_cities.reco.wfm_functions import cwf_from_rwf

    run_dir = Path(root) / str(RUN)
    stem = f"run_{RUN}_0000_ldc1_trg0"
    raw = run_dir / "data" / f"{stem}.waveforms.h5"
    irene = run_dir / "irene" / f"{stem}.v2.3.1.20250218.Kr.irene.h5"
    pmts = DataPMT("demopp", RUN)
    gains = pmts.adc_to_pes.to_numpy()
    calibration = CalibVectors(
        channel_id=pmts.ChannelID.to_numpy(), coeff_blr=abs(pmts.coeff_blr.to_numpy()),
        coeff_c=abs(pmts.coeff_c.to_numpy()), adc_to_pes=gains, adc_to_pes_sipm=np.array([]),
        pmt_active=np.flatnonzero(pmts.Active.to_numpy()).tolist(),
    )
    deconvolution = DeconvParams(n_baseline=14000, thr_trigger=5)
    with tables.open_file(raw, "r") as h5:
        rows = np.flatnonzero(h5.root.Run.events.cols.evt_number[:] == EVENT)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row for event {EVENT} in {raw}")
        corrected = cwf_from_rwf(h5.root.RD.pmtrwf, [int(rows[0])], calibration, deconvolution)[0]
    summed = (corrected / gains[:, None]).sum(axis=0)

    windows = {}
    with tables.open_file(irene, "r") as handle:
        for kind in ("S1", "S2"):
            peaks = handle.get_node(f"/PMAPS/{kind}").read_where(f"event == {EVENT}")
            windows[kind] = []
            for peak in np.unique(peaks["peak"]):
                times = peaks[peaks["peak"] == peak]["time"]
                start, end = float(times.min()) / 1000, float(times.max()) / 1000
                windows[kind].append((start, end if end > start or kind == "S2" else end + SAMPLE_US))
    return np.arange(len(summed)) * SAMPLE_US, summed, windows


def export_waveform(root, path):
    time, charge, windows = calibrate_waveform(root)
    # Preserve just the time range displayed in the preprint.
    samples = (time >= 320) & (time <= 480)
    np.savez_compressed(path, time_us=time[samples], charge_pe=charge[samples],
                        s1_windows_us=windows["S1"], s2_windows_us=windows["S2"])


def plot_waveform(data_path, path):
    from .figures import plt, save

    with np.load(data_path) as data:
        time, charge = data["time_us"], data["charge_pe"]
        windows = {"S1": data["s1_windows_us"], "S2": data["s2_windows_us"]}
    if not len(windows["S1"]) or not len(windows["S2"]):
        raise ValueError(f"Expected S1 and S2 windows for event {EVENT}")
    fig, ax = plt.subplots(figsize=(9, 3), layout="constrained")
    ax.plot(time, charge, color="black", lw=0.8)
    for kind, color in (("S1", "tab:red"), ("S2", "tab:blue")):
        for index, (x0, x1) in enumerate(windows[kind]):
            ax.axvspan(x0, x1, color=color, alpha=0.3, label=kind if index == 0 else None)
    start, end = windows["S1"][0]
    zoom = (start - 2, end + 2)
    samples = (time >= zoom[0]) & (time <= zoom[1])
    inset = ax.inset_axes([0.19, 0.38, 0.26, 0.51])
    inset.plot(time[samples], charge[samples], color="black", lw=0.8)
    inset.axvspan(start, end, color="tab:red", alpha=0.3)
    inset.set(xlim=zoom, ylim=(-1, 1.3 * charge[samples].max()))
    ax.indicate_inset_zoom(inset, edgecolor="0.5")
    ax.set(xlim=(320, 480), xlabel="Time (μs)", ylabel="Amplitude (p.e./bin)")
    ax.legend()
    save(fig, path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export Figure 2 samples in the historical IC environment")
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Path for the .npz export")
    args = parser.parse_args()
    export_waveform(args.analysis_root, args.output)
