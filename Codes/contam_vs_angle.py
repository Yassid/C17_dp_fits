#!/usr/bin/env python3
"""Contaminant (0.45 MeV) fitted yield vs theta_cm, from the contam variant.

Shows that the 0.45 MeV excess is a back-angle feature: its per-bin amplitude
is pinned near the fit floor for theta_cm < 30 deg and carries ~all of its
strength at the back.  Reads results/fine_contam/yields.csv.
"""
import os
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(REPO, "results")
PLOTS = os.path.join(REPO, "plots")

rows = [r for r in csv.DictReader(open(os.path.join(RESULTS, "fine_contam",
                                                    "yields.csv")))
        if r["state"].lower() == "contaminant"]
th = np.array([float(r["bin_center"]) for r in rows])
y = np.array([float(r["counts"]) for r in rows])
e = np.array([float(r["counts_err"]) for r in rows])
# clip absurd (unconstrained, floored) error bars for display
e_disp = np.clip(e, 0, 30)

plt.rcParams.update({"font.size": 13, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True,
                     "ytick.right": True})
fig, ax = plt.subplots(figsize=(8.0, 5.2))
fwd = th < 30
ax.bar(th[fwd], y[fwd], width=2.2, color="#9aa0a6", label=r"$\theta_{\rm c.m.}<30^\circ$ (floored ~0)")
ax.bar(th[~fwd], y[~fwd], width=2.2, color="#6f1d46", label=r"$\theta_{\rm c.m.}\geq30^\circ$")
ax.errorbar(th, y, yerr=e_disp, fmt="none", ecolor="k", capsize=3, lw=1, alpha=0.7)
ax.axvline(30, ls=":", color="k", lw=1.2)
frac = y[~fwd].sum() / y.sum() * 100
ax.set_xlabel(r"$\theta_{\rm c.m.}$ (deg)")
ax.set_ylabel("contaminant fitted yield (counts)")
ax.set_title(r"0.45 MeV contaminant is a back-angle feature "
             f"(${frac:.0f}\\%$ of yield at $\\theta_{{\\rm c.m.}}\\geq30^\\circ$)")
ax.legend(frameon=False)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
p = os.path.join(PLOTS, "contam_vs_angle.png")
fig.savefig(p, dpi=140)
print("wrote", p, f"(back-angle fraction {frac:.1f}%, total {y.sum():.0f} counts)")
