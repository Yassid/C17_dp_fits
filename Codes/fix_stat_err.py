#!/usr/bin/env python3
"""
One-off: correct the per-point statistical error in every dsdo.csv.

The C++ macro previously wrote
    dSigma_stat = sigma_raw * sqrt(1/Y + (dY/Y)^2)
which double-counts the statistics: dY (Minuit GetParError on the fitted
amplitude, from a chi^2 fit to sqrt(N)-weighted histograms) already carries the
Poisson uncertainty of the yield.  The 1/Y term inflates every error bar by up
to sqrt(2).  The macro source is now fixed (dSigma_stat = sigma_raw * dY/Y);
this script brings the already-generated CSVs into line WITHOUT re-running ROOT.

For every results/**/dsdo.csv:
  * read sigma, eff, eff_err from the row and Y, dY from the sibling yields.csv
  * VERIFY: reconstruct the OLD error and require it to match the file value
    (proves correct parsing + keying + formula understanding)
  * rewrite ONLY the dsdo_err_mbsr field with the corrected value
    new_err = sigma * sqrt((dY/Y)^2 + (eff_err/eff)^2)

Rows with Y <= 0 keep dSigma_stat = 0 in the macro, so their error is
unchanged (= sigma*eff_err/eff) and they are left untouched.

Idempotent guard: if a file's errors already match the NEW formula (and not the
old one), it is reported as 'already-fixed' and skipped.
"""
from __future__ import annotations

import os
import csv
import glob
import math

RESULTS = os.path.join(os.path.dirname(__file__), "..", "results")
VERIFY_RTOL = 0.02   # old-formula reconstruction must match file to 2%


def load_yields(path):
    """(bin_lo, bin_hi, state) -> (counts, counts_err)."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r["bin_lo"], r["bin_hi"], r["state"])
            out[key] = (float(r["counts"]), float(r["counts_err"]))
    return out


def old_err(sigma, eff, eff_err, Y, dY):
    if Y > 0:
        return sigma * math.sqrt(1.0 / Y + (dY / Y) ** 2 + (eff_err / eff) ** 2)
    return sigma * (eff_err / eff)


def new_err(sigma, eff, eff_err, Y, dY):
    if Y > 0:
        return sigma * math.sqrt((dY / Y) ** 2 + (eff_err / eff) ** 2)
    return sigma * (eff_err / eff)


def fmt(x):
    return f"{x:.6g}"


def process(dsdo_path):
    d = os.path.dirname(dsdo_path)
    ypath = os.path.join(d, "yields.csv")
    if not os.path.isfile(ypath):
        return ("NO_YIELDS", dsdo_path, 0, 0, 0.0)
    yld = load_yields(ypath)

    with open(dsdo_path) as f:
        lines = f.read().splitlines()
    header = lines[0]
    cols = header.split(",")
    iLo, iHi, iState = cols.index("bin_lo"), cols.index("bin_hi"), cols.index("state")
    iSig = cols.index("dsdo_mbsr")
    iErr = cols.index("dsdo_err_mbsr")
    iEff = cols.index("eff")
    iEffErr = cols.index("eff_err")

    out_lines = [header]
    n_changed = 0
    n_alreadyfixed = 0
    max_resid = 0.0
    worst = None
    sample = None

    for ln in lines[1:]:
        if not ln.strip():
            out_lines.append(ln)
            continue
        tok = ln.split(",")
        key = (tok[iLo], tok[iHi], tok[iState])
        sigma = float(tok[iSig])
        eff = float(tok[iEff])
        eff_err = float(tok[iEffErr])
        cur = float(tok[iErr])
        if key not in yld:
            raise KeyError(f"{dsdo_path}: no yields row for {key}")
        Y, dY = yld[key]

        oe = old_err(sigma, eff, eff_err, Y, dY)
        ne = new_err(sigma, eff, eff_err, Y, dY)

        # idempotency: already fixed?
        if cur > 0 and abs(cur - ne) / max(ne, 1e-30) < VERIFY_RTOL \
                and abs(cur - oe) / max(oe, 1e-30) > VERIFY_RTOL:
            n_alreadyfixed += 1
            out_lines.append(ln)
            continue

        # verify our understanding reproduces the EXISTING (old-formula) value
        if cur > 0:
            resid = abs(cur - oe) / max(oe, 1e-30)
            if resid > max_resid:
                max_resid, worst = resid, (key, cur, oe)
            if resid > VERIFY_RTOL:
                raise AssertionError(
                    f"{dsdo_path} {key}: old-formula reconstruction "
                    f"{oe:.6g} != file {cur:.6g} (resid {resid:.3%}) -- "
                    "aborting, will not rewrite this file")

        if abs(ne - cur) > 0:
            tok[iErr] = fmt(ne)
            n_changed += 1
            if sample is None:
                sample = (key, cur, ne)
        out_lines.append(",".join(tok))

    # write back only if something changed
    if n_changed:
        with open(dsdo_path, "w") as f:
            f.write("\n".join(out_lines) + "\n")

    tag = "FIXED" if n_changed else ("ALREADY" if n_alreadyfixed else "NOCHANGE")
    return (tag, dsdo_path, n_changed, max_resid, sample)


def main():
    paths = sorted(glob.glob(os.path.join(RESULTS, "**", "dsdo.csv"),
                             recursive=True))
    print(f"# {len(paths)} dsdo.csv files found")
    tot_changed = 0
    glob_resid = 0.0
    for p in paths:
        tag, path, nch, resid, sample = process(p)
        tot_changed += nch
        glob_resid = max(glob_resid, resid)
        rel = os.path.relpath(path, RESULTS)
        if tag == "FIXED":
            k, old_v, new_v = sample
            print(f"  {tag:8s} {rel:55s} {nch:3d} rows  "
                  f"e.g. {k[2]}@{k[0]}: {old_v:.4g} -> {new_v:.4g}")
        else:
            print(f"  {tag:8s} {rel:55s}")
    print(f"# total rows corrected: {tot_changed}")
    print(f"# max old-formula reconstruction residual: {glob_resid:.3%} "
          f"(tol {VERIFY_RTOL:.0%}) -- parsing/keying validated")


if __name__ == "__main__":
    main()
