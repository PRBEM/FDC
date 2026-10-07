#!/usr/bin/env python3
"""
Compare two FDC .plt files on the same (E, PA) grid, e.g. a new build
against a reference run.

Usage:  python compare_plt.py new.plt reference.plt
"""
import sys
import numpy as np

NAMES = ["Daa", "Dpp", "Dap"]
THRESHOLDS = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2]


def read_plt(fname):
    """Keep lines with 6 numeric columns (L, E, PA, Daa, Dpp, Dap), sorted by (E, PA)."""
    rows = []
    with open(fname) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 6:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    d = np.array(rows)
    return d[np.lexsort((d[:, 2], d[:, 1]))]


new = read_plt(sys.argv[1])
ref = read_plt(sys.argv[2])
if new.shape != ref.shape or not np.allclose(new[:, 1:3], ref[:, 1:3]):
    sys.exit("ERROR: the two files do not have the same (E, PA) grid")

print(f"{len(new)} grid points\n")
header = ("        cells   zero mism   NaN new/ref  " +
          "  ".join(f">{t:.0e}" for t in THRESHOLDS) + "   max rel diff")
print(header)
for i, name in enumerate(NAMES):
    a, b = new[:, 3 + i], ref[:, 3 + i]
    nan_a, nan_b = np.sum(~np.isfinite(a)), np.sum(~np.isfinite(b))
    finite = np.isfinite(a) & np.isfinite(b)
    zero_mismatch = np.sum(((a == 0) != (b == 0)) & finite)
    # ignore cells that are tiny compared with the grid maximum (near-cancellation
    # in Dap); nanmax so that NaN cells do not wipe out the whole mask
    scale = np.nanmax(np.abs(b[finite])) if finite.any() else 0.0
    mask = finite & (np.abs(b) > 1.0e-6 * scale)
    rel = np.abs(a[mask] - b[mask]) / np.abs(b[mask])
    # after computing rel in compare_plt.py
    print(f"   median {np.median(rel):.1e}, 90th {np.percentile(rel,90):.1e}, 99th {np.percentile(rel,99):.1e}")
    counts = "  ".join(f"{np.sum(rel > t):6d}" for t in THRESHOLDS)
    worst = rel.max() if rel.size else 0.0
    print(f"{name}  {mask.sum():8d}   {zero_mismatch:7d}   {nan_a:5d}/{nan_b:<5d}  {counts}   {worst:.2e}")
    if rel.size and worst > 1e-6:
        j = np.argmax(rel)
        print(f"      worst cell: E = {new[mask][j, 1]:.4g} keV, PA = {new[mask][j, 2]:.2f} deg")