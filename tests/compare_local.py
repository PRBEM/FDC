#!/usr/bin/env python3
"""
Test C: FDC local mode (lam_min = lam_max = 0) against the tiny-window
bounce-average emulation (lam_min = 0, lam_max = LAM_MAX_DEG), converted
to local D via

    D_local = D_emu * cos(a) * (1.30 - 0.56 sin(a)) / lam_max[rad]

Both files are in 1/day, so no 86400 factor appears.

Usage:  python compare_local.py local.plt emulated.plt
"""
import sys
import numpy as np

LAM_MAX_DEG = 0.01     # lam_max used in the emulation run
TOL = 1.0e-5           # pass threshold on max relative difference
NAMES = ["Daa", "Dpp", "Dap"]


def read_plt(fname):
    """Read an FDC .plt file: keep only lines with 6 numeric columns
    (L, E, PA, Daa, Dpp, Dap), sorted by (E, PA)."""
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


loc = read_plt(sys.argv[1])
emu = read_plt(sys.argv[2])

if loc.shape != emu.shape or not np.allclose(loc[:, 1:3], emu[:, 1:3]):
    sys.exit("ERROR: the two files do not have the same (E, PA) grid")

lam = np.radians(LAM_MAX_DEG)
a = np.radians(emu[:, 2])
fac = np.cos(a) * (1.30 - 0.56 * np.sin(a)) / lam

print(f"{len(loc)} grid points\n")
all_ok = True
for i, name in enumerate(NAMES):
    d_loc = loc[:, 3 + i]
    d_exp = emu[:, 3 + i] * fac

    # cells where exactly one of the two is zero: should not happen
    zero_mismatch = np.sum((d_loc == 0) != (d_exp == 0))

    # relative difference, ignoring cells that are tiny compared with the
    # grid maximum (near-cancellation in Dap would inflate the ratio)
    scale = np.max(np.abs(d_exp))
    mask = np.abs(d_exp) > 1.0e-6 * scale
    rel = np.abs(d_loc[mask] - d_exp[mask]) / np.abs(d_exp[mask])
    worst = rel.max() if rel.size else 0.0

    ok = (zero_mismatch == 0) and (worst < TOL)
    all_ok &= ok
    print(f"{name}: max rel diff = {worst:.2e} over {mask.sum()} cells, "
          f"zero/non-zero mismatches = {zero_mismatch}  "
          f"{'PASS' if ok else 'FAIL'}")
    if rel.size and not ok:
        j = np.argmax(rel)
        e, pa = loc[mask][j, 1], loc[mask][j, 2]
        print(f"      worst cell: E = {e:.4g} keV, PA = {pa:.2f} deg")

print("\nTEST C", "PASSED" if all_ok else "FAILED")

# How widespread are the differences?
for i, name in enumerate(NAMES):
    d_loc, d_exp = loc[:, 3+i], emu[:, 3+i] * fac
    m = np.abs(d_exp) > 1e-6 * np.max(np.abs(d_exp))
    rel = np.abs(d_loc[m] - d_exp[m]) / np.abs(d_exp[m])
    print(f"{name}: >1e-5: {(rel > 1e-5).sum()}, >1e-3: {(rel > 1e-3).sum()} of {m.sum()}")

# Dpp: if the problem is geometric, loc/emu depends on PA only, not on E
ratio = np.where(emu[:, 4] != 0, loc[:, 4] / (emu[:, 4] * fac), np.nan)
for pa in np.unique(loc[:, 2])[:10]:
    r = ratio[loc[:, 2] == pa]
    print(f"PA {pa:5.1f}: Dpp loc/emu from {np.nanmin(r):.4f} to {np.nanmax(r):.4f}")
