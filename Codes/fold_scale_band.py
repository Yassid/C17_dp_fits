#!/usr/bin/env python3
"""Fold the absolute-normalization scale uncertainty into the published bands.

The d+16C elastic check (om_elastic.py) gives a cross-potential spread of
~0.86-0.98 in data/OM, i.e. a GLOBAL multiplicative normalization uncertainty
of ~12%.  It is fully correlated across all states and angles (it shifts the
absolute scale, not the relative state-to-state / angle-to-angle picture), so
it must NOT be added per-point in quadrature with the statistical scatter --
it is carried as a separate, correlated outer band.

This writes, for the record:
  results/dsdo_bands_fine.csv  -- per (state, theta): nominal, stat, syst,
                                  point-band = sqrt(stat^2+syst^2), and the
                                  correlated scale band (+-12%).
  results/yield_bands_fine.csv -- per state: nominal yield, syst envelope,
                                  and the same +-12% scale band.
"""
import os
import csv
from collections import defaultdict

REPO = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(REPO, "results")

SCALE = 0.12          # global normalization uncertainty (OMP spread 0.86-0.98)


def main():
    # nominal dsdo + stat error
    dsdo = {}
    with open(os.path.join(RESULTS, "fine", "dsdo.csv")) as f:
        for r in csv.DictReader(f):
            dsdo[(r["state"], round(float(r["bin_center"]), 2))] = (
                float(r["dsdo_mbsr"]), float(r["dsdo_err_mbsr"]))

    # systematic half-range from the 82-variant band
    syst = {}
    with open(os.path.join(RESULTS, "systematics_dsdo_band_fine.csv")) as f:
        for r in csv.DictReader(f):
            syst[(r["state"], round(float(r["bin_center"]), 2))] = \
                float(r["syst_half_range"])

    out = os.path.join(RESULTS, "dsdo_bands_fine.csv")
    n = 0
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "bin_center", "dsdo_mbsr", "stat_err",
                    "syst_half", "scale_err_corr", "band_pointwise",
                    "lo_pointwise", "hi_pointwise", "scale_lo", "scale_hi"])
        for (st, th), (v, stat) in sorted(dsdo.items()):
            sy = syst.get((st, th), 0.0)
            point = (stat**2 + sy**2) ** 0.5          # uncorrelated, per-point
            scl = SCALE * v                            # correlated, global
            w.writerow([st, th, f"{v:.6g}", f"{stat:.6g}", f"{sy:.6g}",
                        f"{scl:.6g}", f"{point:.6g}",
                        f"{max(v-point,0):.6g}", f"{v+point:.6g}",
                        f"{v*(1-SCALE):.6g}", f"{v*(1+SCALE):.6g}"])
            n += 1
    print(f"wrote {out}  ({n} rows; scale = +-{SCALE*100:.0f}% correlated)")

    # yields
    yld_in = os.path.join(RESULTS, "yield_envelope_fine.csv")
    if os.path.isfile(yld_in):
        rows = list(csv.DictReader(open(yld_in)))
        out2 = os.path.join(RESULTS, "yield_bands_fine.csv")
        with open(out2, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["state", "yield", "syst_lo", "syst_hi",
                        "scale_err_corr", "scale_lo", "scale_hi"])
            for r in rows:
                # tolerate either 'nominal' or 'yield' header naming
                nomk = "nominal" if "nominal" in r else list(r)[1]
                v = float(r[nomk])
                lo = r.get("env_lo", r.get("syst_lo", v))
                hi = r.get("env_hi", r.get("syst_hi", v))
                w.writerow([r["state"], f"{v:.6g}", lo, hi,
                            f"{SCALE*v:.6g}", f"{v*(1-SCALE):.6g}",
                            f"{v*(1+SCALE):.6g}"])
        print(f"wrote {out2}  ({len(rows)} states)")
    else:
        print(f"(skipped yields: {yld_in} not found)")


if __name__ == "__main__":
    main()
