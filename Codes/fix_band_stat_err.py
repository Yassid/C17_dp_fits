#!/usr/bin/env python3
"""
Reconcile systematics_dsdo_band_fine.csv with the corrected per-point stat
errors (see fix_stat_err.py).  The band's central values and systematic
envelope (nominal, syst_lo, syst_hi, syst_half_range) depend only on the
*central* dsdo across variants, which did NOT change -- only stat_err did.

So: pull stat_err for each (state, bin_center) straight from the corrected
nominal results/fine/dsdo.csv, and recompute syst_over_stat = syst_half_range /
stat_err.  Everything else is left byte-for-byte.

Sanity assertion: the band's 'nominal' must match the nominal dsdo_mbsr, which
validates the (state, bin_center) keying before any rewrite.
"""
import os
import csv
import math

HERE = os.path.dirname(__file__)
NOMINAL = os.path.join(HERE, "..", "results", "fine", "dsdo.csv")
BAND = os.path.join(HERE, "..", "results", "systematics_dsdo_band_fine.csv")
RTOL = 0.02


def key(state, theta):
    return (state, round(float(theta), 3))


def main():
    # map (state, theta) -> (dsdo, dsdo_err) from corrected nominal
    err = {}
    sig = {}
    with open(NOMINAL) as f:
        for r in csv.DictReader(f):
            k = key(r["state"], r["bin_center"])
            err[k] = float(r["dsdo_err_mbsr"])
            sig[k] = float(r["dsdo_mbsr"])

    with open(BAND) as f:
        rows = list(csv.DictReader(f))
        cols = rows[0].keys() if rows else []

    n_upd = 0
    max_nomresid = 0.0
    for r in rows:
        k = key(r["state"], r["bin_center"])
        if k not in err:
            raise KeyError(f"band row {k} absent from nominal dsdo.csv")
        # validate keying: band 'nominal' must equal nominal dsdo
        nb = float(r["nominal"])
        if sig[k] > 0:
            resid = abs(nb - sig[k]) / sig[k]
            max_nomresid = max(max_nomresid, resid)
            if resid > RTOL:
                raise AssertionError(
                    f"{k}: band nominal {nb:.6g} != dsdo {sig[k]:.6g} "
                    f"({resid:.2%}) -- keying wrong, aborting")
        new_stat = err[k]
        new_ratio = (float(r["syst_half_range"]) / new_stat
                     if new_stat > 0 else 0.0)
        if (abs(new_stat - float(r["stat_err"])) > 0
                or abs(new_ratio - float(r["syst_over_stat"])) > 0):
            n_upd += 1
        r["stat_err"] = f"{new_stat:.6g}"
        r["syst_over_stat"] = f"{new_ratio:.6g}"

    with open(BAND, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cols))
        w.writeheader()
        w.writerows(rows)

    print(f"# band rows updated: {n_upd}/{len(rows)}")
    print(f"# max nominal-keying residual: {max_nomresid:.3%} (tol {RTOL:.0%})")


if __name__ == "__main__":
    main()
