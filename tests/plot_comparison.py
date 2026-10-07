#!/usr/bin/env python3
"""
Comparison figure: FDC local-mode output vs the other code's output.

Top row:    Daa, |Dap|, Dpp [1/s] vs equatorial pitch angle for a few
            energies; FDC as lines, the other code as markers.
Bottom row: ratio FDC / other on the (pitch angle, energy) grid.

Usage:
    python plot_comparison.py fdc.plt L300_i5.dat [--label NAME] [--out NAME]
                              [--energies 0.1 0.3 1 3 10] [--range 0.02]
      --label     legend name for the other code (default: "reference code")
      --out       output file stem (writes .png and .pdf)
      --energies  energies to show in the top row [MeV]
      --range     colour scale of the ratio maps: 1 +/- range (default 0.02)
      --note      extra line for the figure title (e.g. plasma parameters)
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

MC2_CGS = (9.1093897e-28 * 2.99792458e10) ** 2
ME_MEV = 0.51099906


def read_other(fname):
    """Other code: total (n != 0 plus Landau) D/p^2 in 1/s. Columns: E[MeV], PA, Daa, Dap, Dpp."""
    rows, fac, E = [], None, None
    with open(fname) as f:
        for line in f:
            if "omPe=" in line:
                E = float(line.split("E=")[1].split()[0])
                dnorm = float(line.split("=")[-1])
                x = E / ME_MEV
                fac = dnorm / (MC2_CGS * x * (x + 2.0))
                continue
            p = line.split()
            if len(p) != 7:
                continue
            try:
                v = [float(t) for t in p]
            except ValueError:
                continue
            rows.append([E, v[0], (v[1] + v[4]) * fac, (v[2] + v[5]) * fac, (v[3] + v[6]) * fac])
    return np.array(rows)


def read_fdc(fname):
    """FDC .plt: L, E[keV], PA, Daa, Dpp, Dap in 1/day -> E[MeV], PA, Daa, Dap, Dpp in 1/s."""
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
            rows.append([v[1] / 1000.0, v[2], v[3] / 86400.0, v[5] / 86400.0, v[4] / 86400.0])
    return np.array(rows)


def to_grid(d, Es, PAs):
    """Arrays [nE, nPA] for Daa, Dap, Dpp, NaN where a grid point is missing."""
    g = np.full((3, len(Es), len(PAs)), np.nan)
    for e, pa, *vals in d:
        i = np.argmin(np.abs(Es - e) / Es)
        j = np.argmin(np.abs(PAs - pa))
        if abs(Es[i] - e) <= 2e-3 * Es[i] and abs(PAs[j] - pa) < 1e-3:
            g[:, i, j] = vals
    return g


ap = argparse.ArgumentParser()
ap.add_argument("fdc")
ap.add_argument("other")
ap.add_argument("--label", default="JA code")
ap.add_argument("--out", default="fdc_comparison")
ap.add_argument("--energies", type=float, nargs="+", default=[0.1, 0.3, 1.0, 3.0, 10.0])
ap.add_argument("--range", type=float, default=0.02)
ap.add_argument("--note", default="")
args = ap.parse_args()

oth = read_other(args.other)
fdc = read_fdc(args.fdc)
Es = np.unique(np.round(oth[:, 0], 6))
PAs = np.unique(np.round(oth[:, 1], 3))
G_o = to_grid(oth, Es, PAs)
G_f = to_grid(fdc, Es, PAs)

names = [r"$D_{\alpha\alpha}/p^2$", r"$|D_{\alpha p}|/p^2$", r"$D_{pp}/p^2$"]
short = ["Daa", "Dap", "Dpp"]

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
fig, axes = plt.subplots(2, 3, figsize=(14, 8.2), constrained_layout=True,
                         gridspec_kw={"height_ratios": [1.15, 1]})
colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(args.energies)))

# ---- top row: D vs pitch angle ----
for k in range(3):
    ax = axes[0, k]
    for E_req, col in zip(args.energies, colors):
        i = np.argmin(np.abs(np.log(Es / E_req)))
        f = np.abs(G_f[k, i]); o = np.abs(G_o[k, i])
        fm = np.where(f > 0, f, np.nan); om = np.where(o > 0, o, np.nan)
        ax.plot(PAs, fm, "-", color=col, lw=1.6)
        ax.plot(PAs, om, "o", color=col, ms=3.2, mfc="none", mew=0.9)
        # energy label at the right end of the curve
    ax.set_yscale("log")
    top = np.nanmax(np.abs(G_o[k]))
    ax.set_ylim(top * 1e-10, top * 5)          # show ten decades below the maximum
    ax.set_xlim(0, 90)
    ax.set_xlabel("equatorial pitch angle [deg]")
    ax.set_ylabel(names[k] + r"  [s$^{-1}$]")
    ax.set_title(names[k])
    ax.grid(alpha=0.3, which="both", lw=0.4)

# legend: energies (colour) and codes (line style)
from matplotlib.lines import Line2D
h_E = [Line2D([], [], color=c, lw=2) for c in colors]
l_E = [f"{Es[np.argmin(np.abs(np.log(Es / e)))]:.3g} MeV" for e in args.energies]
h_c = [Line2D([], [], color="k", lw=1.6), Line2D([], [], color="k", ls="", marker="o", mfc="none", ms=4)]
l_c = ["FDC", args.label]
axes[0, 0].legend(h_E + h_c, l_E + l_c, fontsize=8.5, loc="upper left", framealpha=0.9)

# ---- bottom row: ratio maps ----
norm = TwoSlopeNorm(vcenter=1.0, vmin=1.0 - args.range, vmax=1.0 + args.range)
dPA = PAs[1] - PAs[0]
PA_edges = np.concatenate([PAs - dPA / 2, [PAs[-1] + dPA / 2]])
lE = np.log10(Es)
dlE = lE[1] - lE[0]
E_edges = 10 ** np.concatenate([lE - dlE / 2, [lE[-1] + dlE / 2]])
for k in range(3):
    ax = axes[1, k]
    f, o = G_f[k], G_o[k]
    ok = (f != 0) & (o != 0) & np.isfinite(f) & np.isfinite(o)
    ratio = np.where(ok, f / np.where(o != 0, o, np.nan), np.nan)
    pm = ax.pcolormesh(PA_edges, E_edges, ratio, cmap="RdBu_r", norm=norm, shading="flat")
    ax.set_facecolor("0.88")
    ax.set_yscale("log")
    ax.set_xlim(0, 90)
    ax.set_xlabel("equatorial pitch angle [deg]")
    ax.set_ylabel("E [MeV]")
    r = ratio[np.isfinite(ratio)]
    within = np.mean(np.abs(r - 1) < 0.01) * 100
    ax.set_title(f"{short[k]}: FDC / {args.label}\n"
                 f"median {np.median(r):.4f}", fontsize=10)
fig.colorbar(pm, ax=axes[1, :], shrink=0.9, label="ratio FDC / " + args.label, extend="both")

title = "Local equatorial diffusion coefficients, L = 3: FDC vs " + args.label
if args.note:
    title += "\n" + args.note
#fig.suptitle(title, fontsize=12)

fig.savefig(args.out + ".png", dpi=160)
fig.savefig(args.out + ".pdf")
print(f"wrote {args.out}.png and {args.out}.pdf")
