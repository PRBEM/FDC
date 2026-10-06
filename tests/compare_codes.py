#!/usr/bin/env python3
"""
Compare diffusion coefficients from several codes on their common grid.

Supported formats (detected from the file extension):
  .plt  FDC local mode: L, E[keV], PA, Daa, Dpp, Dap in 1/day, normalised by p^2
  .dat  Albert's code: per-energy blocks with Dnorm; columns Daan Dapn Dppn Daa0 Dap0 Dpp0
  .h5   QLD++: datasets Daa, Dap_over_p, Dpp_over_p2 in 1/s, plus E_MeV, alpha_deg

Everything is converted to (E [MeV], pitch angle [deg], Daa, Dap/p, Dpp/p^2)
in 1/s. Only grid points present in every file are compared; the first file is
the reference for all ratios.

Usage:
    python compare_codes.py fdc.plt L300_i5.dat ma_local_L3eq.h5 [options]
      --labels FDC Albert QLD++     names used in the output and the figure
      --albert total|n|zero         which columns of the .dat file to use
                                    (total = Daan + Daa0, the default)
      --plot [stem]                 write <stem>.png and <stem>.pdf
      --energies 0.1 1 10           energies for the top row of the figure [MeV]
      --range 0.02                  colour scale of the ratio maps: 1 +/- range
"""
import argparse
import os
import numpy as np

MC2_CGS = (9.1093897e-28 * 2.99792458e10) ** 2
ME_MEV = 0.51099906


def read_fdc(fname, **kw):
    rows = []
    with open(fname) as f:
        for line in f:
            p = line.split()
            if len(p) != 6:
                continue
            try:
                v = [float(t) for t in p]
            except ValueError:
                continue
            # L, E[keV], PA, Daa, Dpp, Dap  ->  E[MeV], PA, Daa, Dap, Dpp  in 1/s
            rows.append([v[1] / 1000.0, v[2], v[3] / 86400.0, v[5] / 86400.0, v[4] / 86400.0])
    return np.array(rows)


def read_albert(fname, albert="total", **kw):
    rows, fac, E = [], None, None
    with open(fname) as f:
        for line in f:
            if "omPe=" in line:
                E = float(line.split("E=")[1].split()[0])
                dnorm = float(line.split("=")[-1])
                x = E / ME_MEV
                fac = dnorm / (MC2_CGS * x * (x + 2.0))       # -> D/p^2 in 1/s
                continue
            p = line.split()
            if len(p) != 7:
                continue
            try:
                v = [float(t) for t in p]
            except ValueError:
                continue
            n_part = np.array(v[1:4])                          # Daan, Dapn, Dppn
            zero_part = np.array(v[4:7])                       # Daa0, Dap0, Dpp0 (Landau)
            d = {"total": n_part + zero_part, "n": n_part, "zero": zero_part}[albert]
            rows.append([E, v[0], d[0] * fac, d[1] * fac, d[2] * fac])
    return np.array(rows)


def read_qldpp(fname, **kw):
    import h5py
    with h5py.File(fname, "r") as f:
        E = f["E_MeV"][:]
        pa = f["alpha_deg"][:]
        daa = f["Daa"][:]                 # [n_alpha, n_E], already 1/s
        dap = f["Dap_over_p"][:]
        dpp = f["Dpp_over_p2"][:]
    rows = []
    for j, e in enumerate(E):
        for i, a in enumerate(pa):
            rows.append([e, a, daa[i, j], dap[i, j], dpp[i, j]])
    return np.array(rows)


READERS = {".plt": read_fdc, ".dat": read_albert, ".h5": read_qldpp, ".hdf5": read_qldpp}


def load(fname, **kw):
    ext = os.path.splitext(fname)[1].lower()
    if ext not in READERS:
        raise SystemExit(f"unknown format: {fname}")
    return READERS[ext](fname, **kw)


def common_grid(datasets):
    """Energies and pitch angles present in every dataset (energies within 0.2%)."""
    PAs = sorted(set.intersection(*[set(np.round(d[:, 1], 3)) for d in datasets]))
    Es = sorted(np.unique(np.round(datasets[0][:, 0], 8)))
    keep = []
    for e in Es:
        if all(np.any(np.abs(d[:, 0] - e) <= 2e-3 * e) for d in datasets[1:]):
            keep.append(e)
    return np.array(keep), np.array(PAs)


def to_grid(d, Es, PAs):
    g = np.full((3, len(Es), len(PAs)), np.nan)
    for e, pa, *vals in d:
        i = int(np.argmin(np.abs(Es - e) / Es))
        j = int(np.argmin(np.abs(PAs - pa)))
        if abs(Es[i] - e) <= 2e-3 * Es[i] and abs(PAs[j] - pa) < 1e-3:
            g[:, i, j] = vals
    return g


ap = argparse.ArgumentParser()
ap.add_argument("files", nargs="+")
ap.add_argument("--labels", nargs="+")
ap.add_argument("--albert", default="total", choices=["total", "n", "zero"])
ap.add_argument("--plot", nargs="?", const="code_comparison", default=None)
ap.add_argument("--energies", type=float, nargs="+", default=[0.1, 1.0, 10.0])
ap.add_argument("--range", type=float, default=0.02)
ap.add_argument("--note", default="")
args = ap.parse_args()

labels = args.labels or [os.path.basename(f) for f in args.files]
if len(labels) != len(args.files):
    raise SystemExit("--labels must give one name per file")

data = [load(f, albert=args.albert) for f in args.files]
Es, PAs = common_grid(data)
grids = [to_grid(d, Es, PAs) for d in data]
print(f"common grid: {len(Es)} energies ({Es.min():.4g} - {Es.max():.4g} MeV), "
      f"{len(PAs)} pitch angles\n")

names = ["Daa", "Dap/p", "Dpp/p2"]
ref, ref_label = grids[0], labels[0]
for k, name in enumerate(names):
    print(f"--- {name}   (ratio to {ref_label})")
    print("                        both>0   only this  only ref    median     "
          "5%-95%            |r-1|<10%   <1%")
    for g, lab in zip(grids[1:], labels[1:]):
        a, b = g[k], ref[k]
        scale = np.nanmax(np.abs(b))
        sig_a = np.abs(a) > 1e-6 * scale
        sig_b = np.abs(b) > 1e-6 * scale
        both = sig_a & sig_b & np.isfinite(a) & np.isfinite(b)
        r = a[both] / b[both]
        if r.size:
            p5, p50, p95 = np.percentile(r, [5, 50, 95])
            print(f"  {lab:<20} {both.sum():8d} {np.sum(sig_a & ~sig_b):10d} "
                  f"{np.sum(~sig_a & sig_b):10d}   {p50:8.4f}   {p5:7.3f} - {p95:7.3f}"
                  f"      {np.mean(np.abs(r-1) < 0.10)*100:6.1f}%  {np.mean(np.abs(r-1) < 0.01)*100:6.1f}%")
        else:
            print(f"  {lab:<20} no overlapping cells")
    print()

if args.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.lines import Line2D

    nrat = len(grids) - 1
    fig, axes = plt.subplots(1 + nrat, 3, figsize=(14, 4.1 * (1 + nrat)),
                             constrained_layout=True, squeeze=False)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(args.energies)))
    styles = ["-", "o", "s", "^"]

    for k in range(3):
        ax = axes[0, k]
        for E_req, col in zip(args.energies, colors):
            i = int(np.argmin(np.abs(np.log(Es / E_req))))
            for m, (g, lab) in enumerate(zip(grids, labels)):
                y = np.abs(g[k, i])
                y = np.where(y > 0, y, np.nan)
                if m == 0:
                    ax.plot(PAs, y, "-", color=col, lw=1.6)
                else:
                    ax.plot(PAs, y, styles[m % len(styles)], color=col, ms=3.5,
                            mfc="none", mew=0.9)
        ax.set_yscale("log")
        top = np.nanmax(np.abs(ref[k]))
        ax.set_ylim(top * 1e-10, top * 5)
        ax.set_xlim(0, 90)
        ax.set_xlabel("pitch angle [deg]")
        ax.set_ylabel(names[k] + r"  [s$^{-1}$]")
        ax.set_title(names[k])
        ax.grid(alpha=0.3, which="both", lw=0.4)
    h = [Line2D([], [], color=c, lw=2) for c in colors] + \
        [Line2D([], [], color="k", lw=1.6)] + \
        [Line2D([], [], color="k", ls="", marker=styles[m % len(styles)], mfc="none", ms=4)
         for m in range(1, len(labels))]
    l = [f"{Es[int(np.argmin(np.abs(np.log(Es/e))))]:.3g} MeV" for e in args.energies] + labels
    axes[0, 0].legend(h, l, fontsize=8, loc="lower right", framealpha=0.9)

    norm = TwoSlopeNorm(vcenter=1.0, vmin=1.0 - args.range, vmax=1.0 + args.range)
    dPA = PAs[1] - PAs[0]
    PA_edges = np.concatenate([PAs - dPA / 2, [PAs[-1] + dPA / 2]])
    lE = np.log10(Es)
    dlE = np.median(np.diff(lE))
    E_edges = 10 ** np.concatenate([lE - dlE / 2, [lE[-1] + dlE / 2]])
    for m in range(nrat):
        for k in range(3):
            ax = axes[1 + m, k]
            a, b = grids[m + 1][k], ref[k]
            ok = (a != 0) & (b != 0) & np.isfinite(a) & np.isfinite(b)
            ratio = np.where(ok, a / np.where(b != 0, b, np.nan), np.nan)
            pm = ax.pcolormesh(PA_edges, E_edges, ratio, cmap="RdBu_r", norm=norm,
                               shading="flat")
            ax.set_facecolor("0.88")
            ax.set_yscale("log")
            ax.set_xlim(0, 90)
            ax.set_xlabel("pitch angle [deg]")
            ax.set_ylabel("E [MeV]")
            r = ratio[np.isfinite(ratio)]
            ax.set_title(f"{names[k]}: {labels[m+1]} / {ref_label}\n"
                         f"median {np.median(r):.4f},  "
                         f"{np.mean(np.abs(r-1) < 0.01)*100:.1f}% of cells within 1%",
                         fontsize=10)
        fig.colorbar(pm, ax=axes[1 + m, :], shrink=0.9, extend="both",
                     label=f"ratio / {ref_label}")

    if args.note:
        fig.suptitle(args.note, fontsize=12)
    fig.savefig(args.plot + ".png", dpi=150)
    fig.savefig(args.plot + ".pdf")
    print(f"wrote {args.plot}.png and {args.plot}.pdf")
