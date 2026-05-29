#!/usr/bin/env python3
"""Inclusive-fit residuals from results/<variant>/inclusive_fit.csv.

Top panel: data and fit components (stacked colors), zoomed 0-3 MeV.
Mid panel: residual = (data - fit_total) / data_err per bin.
Bottom row: integrated excess (data - fit_total) in 0.7-2.5 MeV.

Default variant = "fine" (nominal). Override with --variant.
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]


def load(variant: str) -> dict[str, np.ndarray]:
    path = REPO / "results" / variant / "inclusive_fit.csv"
    if not path.exists():
        raise SystemExit(f"missing: {path}  (run C16_pd_AngBins.C first)")
    rows = list(csv.DictReader(path.open()))
    keys = rows[0].keys()
    return {k: np.array([float(r[k]) for r in rows]) for k in keys}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="fine")
    ap.add_argument("--xmax", type=float, default=3.0)
    ap.add_argument("--contam-lo", type=float, default=0.7)
    ap.add_argument("--contam-hi", type=float, default=2.5)
    args = ap.parse_args()

    d = load(args.variant)
    Ex = d["Ex_MeV"]
    data, derr = d["data"], d["data_err"]
    tot = d["fit_total"]
    gauss, bw = d["sum_gauss"], d["sum_bw"]
    ps1n, ps2n, bg = d["ps1n"], d["ps2n"], d["bg"]

    resid = np.where(derr > 0, (data - tot) / derr, 0.0)
    # safe symmetric pull: skip empty bins (derr=0 typical at <-0.5 MeV)

    fig, (ax_top, ax_mid) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.05),
    )

    # Top: data + components
    sel = Ex <= args.xmax + 0.5
    ax_top.errorbar(Ex[sel], data[sel], yerr=derr[sel], fmt="o",
                    color="k", ms=3, capsize=0, label="Data")
    ax_top.plot(Ex[sel], tot[sel], "r-", lw=2, label="Total fit")
    ax_top.plot(Ex[sel], gauss[sel], color="C0", lw=1, label="Bound (3 Gauss)")
    ax_top.plot(Ex[sel], bw[sel], color="0.3", lw=1, label="Unbound BWs")
    ax_top.plot(Ex[sel], ps1n[sel], color="C2", lw=1.2, label="1n phase space")
    ax_top.plot(Ex[sel], ps2n[sel], color="C2", lw=1, ls="--", label="2n phase space")
    ax_top.plot(Ex[sel], bg[sel], color="0.6", lw=1, ls=":", label="Linear bg")
    ax_top.set_xlim(-0.5, args.xmax)
    ax_top.set_ylabel("Counts")
    ax_top.set_title(f"Inclusive fit — {args.variant}")
    ax_top.legend(loc="upper left", fontsize=8, ncols=2)
    ax_top.grid(True, alpha=0.3)

    # Mid: residual
    sel_r = sel & (derr > 0)
    ax_mid.axhline(0, color="k", lw=0.8)
    ax_mid.axhspan(-1, 1, color="0.85", alpha=0.5)
    ax_mid.errorbar(Ex[sel_r], resid[sel_r], yerr=1.0, fmt="o",
                    color="k", ms=3, capsize=0)
    # mark the contaminant window
    ax_mid.axvspan(args.contam_lo, args.contam_hi, color="orange", alpha=0.15)
    ax_mid.set_xlabel("Excitation energy (MeV)")
    ax_mid.set_ylabel("(data $-$ fit) / err")
    ax_mid.grid(True, alpha=0.3)
    ax_mid.set_ylim(-5, 5)

    # Integrated excess in the contaminant window
    mask = (Ex >= args.contam_lo) & (Ex <= args.contam_hi)
    excess = (data - tot)[mask].sum()
    excess_err = np.sqrt((derr[mask]**2).sum())
    n_data = data[mask].sum()
    n_model = tot[mask].sum()

    txt = (f"Window {args.contam_lo:.2f}–{args.contam_hi:.2f} MeV\n"
           f"  $\\sum$ data  = {n_data:.0f}\n"
           f"  $\\sum$ model = {n_model:.0f}\n"
           f"  excess = {excess:+.0f} $\\pm$ {excess_err:.0f}  "
           f"({excess/excess_err:+.1f} $\\sigma$)")
    ax_top.text(0.98, 0.05, txt, transform=ax_top.transAxes,
                ha="right", va="bottom", fontsize=8,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.6"))

    out = REPO / "plots" / f"inclusive_residuals_{args.variant}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"window [{args.contam_lo}, {args.contam_hi}] MeV:")
    print(f"  data  = {n_data:.0f}")
    print(f"  model = {n_model:.0f}")
    print(f"  excess = {excess:+.0f} +/- {excess_err:.0f}  "
          f"({excess/excess_err:+.1f} sigma)")


if __name__ == "__main__":
    main()
