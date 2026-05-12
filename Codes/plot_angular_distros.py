"""
Read dsdo.csv (produced by C16_pd_AngBins.C) and plot dsigma/dOmega vs theta_CM
for each of the 10 17C states in a 2x5 grid.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "results", "dsdo.csv"))

# Preserve the natural state order (rows in yields.csv keep it).
state_order = list(df["state"].drop_duplicates())

state_titles = {
    "gs_3half+":       r"g.s. 3/2$^+$",
    "ex_1half+_217":   r"0.217 1/2$^+$  (1d 2s)",
    "ex_5half+_335":   r"0.335 5/2$^+$",
    "ex2.763_1half-":  r"2.763 1/2$^-$",
    "ex2.980_3half+":  r"2.980 3/2$^+$",
    "ex3.661_NEW":     r"3.661 NEW (3/2$^+$/5/2$^+$)",
    "ex4.231_3half+":  r"4.231 3/2$^+$",
    "ex4.841_1half+":  r"4.841 1/2$^+$",
    "ex5.91_3half+":   r"5.91 3/2$^+$",
    "ex6.30_5half+":   r"6.30 5/2$^+$",
}

fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True)

for i, state in enumerate(state_order):
    ax = axes[i // 5, i % 5]
    sub = df[df["state"] == state].sort_values("bin_center")
    ax.errorbar(sub["bin_center"], sub["dsdo_mbsr"], yerr=sub["dsdo_err_mbsr"],
                fmt="o-", color="C0", ms=7, capsize=3)
    ax.set_yscale("log")
    ax.set_xlim(8, 42)
    ax.set_xticks([10, 15, 20, 25, 30, 35, 40])
    ax.set_ylim(1e-3, 20)
    ax.set_title(state_titles.get(state, state), fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    if i % 5 == 0:
        ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")
    if i // 5 == 1:
        ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")

fig.suptitle("16C(d,p)17C — angular distributions extracted from per-θ_CM Ex spectral fits",
             fontsize=13)
fig.tight_layout()
out = os.path.join(HERE, "..", "plots", "plots_angular_distros.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"saved: {out}")

# Print a summary table
print("\nSummary (dsigma/dOmega in mb/sr):")
piv = df.pivot_table(index="bin_center", columns="state", values="dsdo_mbsr")
piv = piv[state_order]
print(piv.round(3).to_string())
