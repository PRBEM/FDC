"""
Shared reader for the int_N_omega diagnostic dumps (intn_*.txt).

Handles both the 13-column (pre-join-key) and 15-column (with iroot and
theta_res) formats, recovers the record index when the I3 field overflows
at nint_N_res_cone > 999, and reconstructs N(omega) by the same trapezoidal
rule the Fortran uses.
"""

import math
import collections

COLS15 = ("call i st path nres iroot alpha lam thres theta thmax "
          "om n2 dwdk Ntmp").split()
COLS13 = ("call i st path nres alpha lam theta thmax "
          "om n2 dwdk Ntmp").split()
INTCOLS = ("call", "i", "st", "path", "nres", "iroot")


def load(path, thmax=None):
    """Read one intn file into a list of dicts.

    thmax : if given, keep only records with this theta_max (guards against
            files that had an older run appended to them).
    """
    rows = []
    prev_call, idx = None, 0
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) == 15:
                cols = COLS15
            elif len(t) == 13:
                cols = COLS13
            else:
                continue
            overflow = "***" in t[1]
            if overflow:
                t[1] = "0"
            if any("*" in x for x in t):
                continue
            try:
                d = dict(zip(cols, map(float, t)))
            except ValueError:
                continue
            for k in INTCOLS:
                if k in d:
                    d[k] = int(d[k])
            d.setdefault("iroot", 1)
            d.setdefault("thres", -1.0)
            # rebuild i by position when the field overflowed
            if d["call"] != prev_call:
                prev_call, idx = d["call"], 0
            idx += 1
            if overflow:
                d["i"] = idx
            rows.append(d)
    if thmax is not None:
        rows = [r for r in rows if abs(r["thmax"] - thmax) < 1e-6]
    return rows


def group(rows):
    """Group records into one list per int_N_omega invocation, ordered by i.

    Keeps only the first record for each i within a call: files written with
    position='append' can contain a whole earlier run duplicated, which would
    otherwise give a zero step and hence N = 0.
    """
    g = collections.defaultdict(dict)
    for r in rows:
        g[r["call"]].setdefault(r["i"], r)
    return {k: [v[i] for i in sorted(v)] for k, v in g.items()}


def resonance_key(cl):
    """Identity of the resonance a call belongs to, for matching across runs."""
    r = cl[0]
    return (round(r["alpha"], 2), round(r["lam"], 4),
            round(r["thres"], 3), r["nres"], r["iroot"])


def N_of_omega(cl):
    """Trapezoidal reconstruction of N(omega) from the logged integrand."""
    if len(cl) < 2:
        return 0.0
    step = (cl[1]["theta"] - cl[0]["theta"]) * math.pi / 180.0
    tot = 0.0
    for i, r in enumerate(cl):
        w = 0.5 if (i == 0 or i == len(cl) - 1) else 1.0
        tot += w * r["Ntmp"]
    return tot * step


def load_by_resonance(path, thmax=None):
    """{resonance_key: [records]} for one file."""
    out = {}
    for cl in group(load(path, thmax)).values():
        out[resonance_key(cl)] = cl
    return out
