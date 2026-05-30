#!/usr/bin/env python3
"""How much do the NON-RESONANT components affect the extracted dsigma/dOmega,
as a function of angle, and which component is it?

The spectral model has three non-resonant terms, all FREE in every angular bin
(C16_pd_AngBins.C: p[38]=linear bg, p[37]=1n phase space, p[39]=2n phase space;
per-bin fit frees amplitudes + bg, comment "Lock everything ... except
amplitudes + bg").  The systematics campaign scales each by 0.3x / 3x:
  - linear bg : fine_bg_x0.30   / fine_bg_x3.00
  - 1n PS     : fine_ps1n_x0.30 / fine_ps1n_x3.00
  - 2n PS     : fine_ps2n_x0.30 / fine_ps2n_x3.00
The induced fractional change in each unbound state's dsigma/dOmega, binned in
theta_cm, is a direct measure of that component's leverage vs angle.

Outputs bg_vs_angle.csv and bg_vs_angle.png.
"""
import os
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(REPO, "results")
PLOTS = os.path.join(REPO, "plots")

UNBOUND = ["ex2.763_1half-", "ex2.980_3half+", "ex3.661_NEW", "ex4.231_3half+",
           "ex4.841_1half+", "ex5.91_3half+", "ex6.30_5half+"]

COMPONENTS = [
    ("linear bg", "fine_bg_x0.30",   "fine_bg_x3.00",   "#6f1d46", "o-"),
    ("1n phase sp.", "fine_ps1n_x0.30", "fine_ps1n_x3.00", "#2a7fb8", "s-"),
    ("2n phase sp.", "fine_ps2n_x0.30", "fine_ps2n_x3.00", "#2ca02c", "^-"),
]


def load(variant):
    out = {}
    fp = os.path.join(RESULTS, variant, "dsdo.csv")
    with open(fp) as f:
        for r in csv.DictReader(f):
            out[(r["state"], round(float(r["bin_center"]), 2))] = \
                float(r["dsdo_mbsr"])
    return out


def leverage_by_theta(nom, lo, hi):
    by = defaultdict(list)
    for (st, th), v in nom.items():
        if st not in UNBOUND or v <= 0:
            continue
        if (st, th) in lo and (st, th) in hi:
            by[th].append(abs(hi[(st, th)] - lo[(st, th)]) / v)
    th = sorted(by)
    med = [float(np.median(by[t])) for t in th]
    return np.array(th), np.array(med)


def main():
    nom = load("fine")
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    table = {}
    for label, vlo, vhi, col, mk in COMPONENTS:
        try:
            th, med = leverage_by_theta(nom, load(vlo), load(vhi))
        except FileNotFoundError as e:
            print("skip", label, e); continue
        ax.plot(th, 100 * med, mk, color=col, lw=2, ms=6, label=label)
        table[label] = dict(zip(th, med))
        print(f"{label:14s} median leverage by theta:",
              " ".join(f"{t:.0f}:{m*100:.1f}%" for t, m in zip(th, med)))

    ax.axvline(30, ls=":", color="k", lw=1.2)
    ax.text(30.4, ax.get_ylim()[1] * 0.92, r"L-scan cut ($\theta\leq30^\circ$)",
            rotation=90, va="top", fontsize=10)
    ax.set_xlabel(r"$\theta_{\rm c.m.}$ (deg)")
    ax.set_ylabel("median leverage on d$\\sigma$/d$\\Omega$ (%)\n"
                  "(full 0.3$\\times$–3$\\times$ swing / nominal)")
    ax.set_title("Non-resonant components: leverage on the extracted "
                 "cross section vs angle")
    ax.legend(frameon=False, title="component scaled 0.3$\\times$–3$\\times$")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(PLOTS, "bg_vs_angle.png")
    fig.savefig(p, dpi=140)
    print("\nwrote", p)

    # csv
    allth = sorted({t for d in table.values() for t in d})
    with open(os.path.join(RESULTS, "bg_vs_angle.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta_cm"] + [c[0] for c in COMPONENTS])
        for t in allth:
            w.writerow([t] + [round(table.get(c[0], {}).get(t, float("nan")), 4)
                              for c in COMPONENTS])
    print("wrote", os.path.join(RESULTS, "bg_vs_angle.csv"))


if __name__ == "__main__":
    main()
