#!/usr/bin/env python3
"""Figures for the analysis-overview deck.

  1. bound_states_ad.png  -- dsigma/dOmega of the 3 bound states overlaid
  2. unbound_ad_grid.png  -- 2x4 grid of the 7 unbound-state ADs with the
                             82-variant systematic band + corrected stat bars
  3. lscan_bars.png       -- chi2/nu of each L for the 7 unbound states (bars),
                             best-L highlighted

Reads results/dsdo.csv (corrected) and the variant dirs for the band, plus
results/fine/L_scan.csv for the L summary.  Pure matplotlib, repo style.
"""
import os
import csv
import glob
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(REPO, "results")
PLOTS = os.path.join(REPO, "plots")

plt.rcParams.update({
    "font.size": 13, "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.dpi": 140,
})

BOUND = [
    ("gs_3half+",      r"g.s. $3/2^+$",     "#1f77b4", "o"),
    ("ex_1half+_217",  r"0.217 $1/2^+$",    "#d62728", "s"),
    ("ex_5half+_335",  r"0.335 $5/2^+$",    "#2ca02c", "^"),
]

UNBOUND = [
    ("ex2.763_1half-", r"2.763 $1/2^-$"),
    ("ex2.980_3half+", r"2.980 $3/2^+$"),
    ("ex3.661_NEW",    r"3.661 NEW"),
    ("ex4.231_3half+", r"4.231 $3/2^+$"),
    ("ex4.841_1half+", r"4.841 $1/2^+$"),
    ("ex5.91_3half+",  r"5.91 $3/2^+$"),
    ("ex6.30_5half+",  r"6.30 $5/2^+$"),
]


def read_dsdo(path):
    """state -> (theta[], dsdo[], err[]) sorted by theta."""
    out = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                out[r["state"]].append((float(r["bin_center"]),
                                        float(r["dsdo_mbsr"]),
                                        float(r["dsdo_err_mbsr"])))
            except (KeyError, ValueError):
                continue
    res = {}
    for s, v in out.items():
        v.sort()
        a = np.array(v)
        res[s] = (a[:, 0], a[:, 1], a[:, 2])
    return res


def fig_bound():
    d = read_dsdo(os.path.join(RESULTS, "dsdo.csv"))
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    for key, lab, col, mk in BOUND:
        if key not in d:
            continue
        th, ds, er = d[key]
        m = ds > 0
        ax.errorbar(th[m], ds[m], yerr=er[m], fmt=mk, ms=6, color=col,
                    capsize=3, lw=1.4, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta_{\rm c.m.}$ (deg)")
    ax.set_ylabel(r"d$\sigma$/d$\Omega$ (mb/sr)")
    ax.set_title(r"$^{16}$C(d,p)$^{17}$C -- bound-state angular distributions")
    ax.legend(frameon=False, fontsize=12)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    p = os.path.join(PLOTS, "bound_states_ad.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


def fig_unbound_grid():
    base = os.path.join(RESULTS, "fine")
    nominal = read_dsdo(os.path.join(base, "dsdo.csv"))
    # systematic band from all fine* variant dirs (exclude the same set the
    # repo excludes)
    exclude = ("no_3661", "contam_plus_relax")
    variant_dsdo = []
    for d in sorted(glob.glob(os.path.join(RESULTS, "fine*"))):
        if not os.path.isdir(d):
            continue
        tag = os.path.basename(d)
        if any(x in tag for x in exclude):
            continue
        fp = os.path.join(d, "dsdo.csv")
        if os.path.isfile(fp):
            variant_dsdo.append(read_dsdo(fp))

    fig, axes = plt.subplots(2, 4, figsize=(16, 7.6), sharex=True)
    axes = axes.flatten()
    for ax, (key, lab) in zip(axes, UNBOUND):
        if key not in nominal:
            ax.set_visible(False); continue
        th, ds, er = nominal[key]
        lo = np.full_like(ds, np.inf); hi = np.full_like(ds, -np.inf)
        for vd in variant_dsdo:
            if key not in vd:
                continue
            vt, vs, _ = vd[key]
            for i, t in enumerate(th):
                j = np.where(np.isclose(vt, t, atol=0.05))[0]
                if len(j):
                    lo[i] = min(lo[i], vs[j[0]]); hi[i] = max(hi[i], vs[j[0]])
        bad = ~np.isfinite(lo); lo[bad] = ds[bad]; hi[bad] = ds[bad]
        ax.fill_between(th, lo, hi, color="0.82", zorder=1,
                        label="syst. envelope")
        ax.errorbar(th, ds, yerr=er, fmt="o", ms=4.5, color="k", capsize=2,
                    zorder=3, label=r"data $\pm$ stat")
        ax.set_yscale("log")
        ax.set_title(lab, fontsize=12)
        ax.grid(alpha=0.25, which="both")
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    for ax in axes[4:]:
        ax.set_xlabel(r"$\theta_{\rm c.m.}$ (deg)")
    for ax in (axes[0], axes[4]):
        ax.set_ylabel(r"d$\sigma$/d$\Omega$ (mb/sr)")
    if len(UNBOUND) < len(axes):
        for ax in axes[len(UNBOUND):]:
            ax.set_visible(False)
    fig.suptitle(r"$^{16}$C(d,p)$^{17}$C unbound resonances -- "
                 r"d$\sigma$/d$\Omega$ with 82-variant systematic band",
                 y=1.0, fontsize=14)
    fig.tight_layout()
    p = os.path.join(PLOTS, "unbound_ad_grid.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    print("wrote", p)


def fig_lscan():
    rows = list(csv.DictReader(open(os.path.join(RESULTS, "fine", "L_scan.csv"))))
    label = {r[0]: r[1] for r in
             [(k, v) for k, v in zip([u[0] for u in UNBOUND],
                                     [u[1] for u in UNBOUND])]}
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    n = len(rows); w = 0.2
    x = np.arange(n)
    colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52"]
    for li, c in zip(range(4), colors):
        vals = [float(r[f"c2v_L{li}"]) for r in rows]
        ax.bar(x + (li - 1.5) * w, vals, w, color=c, label=f"L={li}")
    # mark best L
    for i, r in enumerate(rows):
        bl = int(r["bestL"])
        v = float(r[f"c2v_L{bl}"])
        ax.plot(i + (bl - 1.5) * w, v, "k*", ms=12, zorder=5)
    ax.axhline(1.0, ls=":", color="0.4", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([label.get(r["state"], r["state"]) for r in rows],
                       rotation=25, ha="right", fontsize=10)
    ax.set_ylabel(r"$\chi^2/\nu$ (AD vs FRESCO L)")
    ax.set_title(r"L assignment: $\chi^2/\nu$ per transferred L "
                 r"($\bigstar$ = best, corrected errors)")
    ax.legend(frameon=False, ncol=4, fontsize=11)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    p = os.path.join(PLOTS, "lscan_bars.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    fig_bound()
    fig_unbound_grid()
    fig_lscan()
