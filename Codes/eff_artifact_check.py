"""
Diagnostic: is the L=0 best-fit for ex5.91 / ex6.30 a genuine forward peak
or an artefact of the Ex-dependent efficiency correction?

Two tests:
  (A) Within-bin efficiency gradient.  The 13.75 deg analysis bin spans
      12.5-15.0 deg.  If eff varies by a large factor across that width,
      a single bin-center eff is meaningless and the corrected point is
      unreliable.
  (B) L-scan vs eff-floor.  Re-run the single-norm L=0..3 fit while
      raising the eff-floor.  If L=0 only wins when the low-eff forward
      bins are included, the flip is correction-driven.

Run with /home/yassid/gnn_env/bin/python.
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
FRESCO = os.path.join(HERE, "..", "fresco", "outputs")
EFFDIR = "/home/yassid/C16_dp/C16_dp/efficiency/"
RCSDIR = "/home/yassid/fair_install/16C_dp/RCS/"

# ---- load the 8 efficiency tables (same set as C16_pd_AngBins.C) -----------
eff_srcs = [
    (0.0, RCSDIR + "efficiency_16Cdp_0m.txt"),
    (2.0, EFFDIR + "efficiency_16Cdp_2m.txt"),
    (3.0, EFFDIR + "efficiency_16Cdp_3m.txt"),
    (4.0, EFFDIR + "efficiency_16Cdp_4m.txt"),
    (4.8, EFFDIR + "efficiency_16Cdp_4.8m.txt"),
    (5.0, EFFDIR + "efficiency_16Cdp_5m.txt"),
    (5.9, EFFDIR + "efficiency_16Cdp_5.9m.txt"),
    (6.0, EFFDIR + "efficiency_16Cdp_6m.txt"),
]
eff_tables = []
for ex, path in eff_srcs:
    a, e, de = np.loadtxt(path, comments="#").T
    eff_tables.append((ex, interp1d(a, e, bounds_error=False,
                                    fill_value=(e[0], e[-1]))))
eff_tables.sort(key=lambda t: t[0])
eff_Ex = np.array([t[0] for t in eff_tables])


def eff2d(Ex, theta):
    """Bilinear (Ex, theta) interpolation, clamped to the table Ex range."""
    if Ex <= eff_Ex[0]:
        return float(eff_tables[0][1](theta))
    if Ex >= eff_Ex[-1]:
        return float(eff_tables[-1][1](theta))
    k = np.searchsorted(eff_Ex, Ex) - 1
    g = (Ex - eff_Ex[k]) / (eff_Ex[k + 1] - eff_Ex[k])
    return float((1 - g) * eff_tables[k][1](theta) + g * eff_tables[k + 1][1](theta))


# ---- TEST A: within-bin efficiency gradient -------------------------------
print("=" * 72)
print("TEST A -- efficiency across the 13.75 deg analysis bin (12.5-15.0 deg)")
print("=" * 72)
for label, Ex in [("ex5.91 (Ex=5.98)", 5.975), ("ex6.30 (Ex=6.36)", 6.355)]:
    e_lo = eff2d(Ex, 12.5)
    e_ct = eff2d(Ex, 13.75)
    e_hi = eff2d(Ex, 15.0)
    print(f"  {label}:  eff(12.5)={e_lo:.3f}  eff(13.75)={e_ct:.3f}  "
          f"eff(15.0)={e_hi:.3f}   spread x{e_hi / e_lo:.1f}")
print("  -> a single bin-center eff cannot represent a x3-4 swing; the")
print("     13.75 deg corrected point is dominated by interpolation choice.\n")

# ---- TEST B: L-scan vs eff-floor ------------------------------------------
df = pd.read_csv(os.path.join(HERE, "..", "results", "dsdo.csv"))
curves = {}
for L in range(4):
    for s in range(1, 8):
        th, sig = np.loadtxt(os.path.join(FRESCO, f"dp17C_L{L}_state{s}.dat")).T
        curves[(L, s)] = interp1d(th, sig, bounds_error=False, fill_value=np.nan)

KMAX = 30.0
SIGFLOOR = 0.03


def fit_one(theta, sig, err, curve):
    m = curve(theta)
    W = 1.0 / err ** 2
    den = float(np.nansum(m * m * W))
    if den <= 0:
        return 0.0, float("nan")
    a = float(np.nansum(sig * m * W)) / den
    chi2 = float(np.nansum(((sig - a * m) / err) ** 2))
    ndf = max(len(theta) - 1, 1)
    return a, chi2 / ndf


# state name -> FRESCO state index (1-based)
states = [("ex5.91_3half+", 6), ("ex6.30_5half+", 7), ("ex4.841_1half+", 5)]
L_lbl = ["L0", "L1", "L2", "L3"]

for name, sidx in states:
    print("=" * 72)
    print(f"TEST B -- {name}  (FRESCO state {sidx})")
    print("=" * 72)
    sub = df[df["state"] == name].sort_values("bin_center")
    print(f"{'eff-floor':>9s} | {'N':>2s} | "
          + " | ".join(f"{l:>10s}" for l in L_lbl) + " | best")
    print("-" * 72)
    for floor in [0.15, 0.30, 0.50, 0.70]:
        m = ((sub["dsdo_mbsr"] > SIGFLOOR) & (sub["dsdo_err_mbsr"] > 0)
             & (sub["bin_center"] <= KMAX) & (sub["eff"] >= floor))
        t = sub["bin_center"].to_numpy()[m]
        s = sub["dsdo_mbsr"].to_numpy()[m]
        e = sub["dsdo_err_mbsr"].to_numpy()[m]
        if len(t) < 2:
            print(f"{floor:9.2f} | {len(t):2d} | (too few points)")
            continue
        cells, c2s = [], []
        for L in range(4):
            a, c2v = fit_one(t, s, e, curves[(L, sidx)])
            c2s.append(c2v if (np.isfinite(c2v) and a > 0) else np.inf)
            cells.append(f"{c2v:10.2f}" if np.isfinite(c2v) else f"{'nan':>10s}")
        best = int(np.argmin(c2s))
        print(f"{floor:9.2f} | {len(t):2d} | " + " | ".join(cells)
              + f" | L={best}")
    print()

# ---- summary plot ---------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
theta_grid = np.linspace(10, 32, 400)
for ax, (name, sidx, title) in zip(
        axes, [("ex5.91_3half+", 6, r"5.91 3/2$^+$"),
               ("ex6.30_5half+", 7, r"6.30 5/2$^+$")]):
    sub = df[df["state"] == name].sort_values("bin_center")
    t = sub["bin_center"].to_numpy()
    s = sub["dsdo_mbsr"].to_numpy()
    e = sub["dsdo_err_mbsr"].to_numpy()
    ef = sub["eff"].to_numpy()
    bad = np.isclose(t, 13.75)                       # the unreliable bin
    good = (~bad) & (s > SIGFLOOR) & (e > 0) & (t <= KMAX)
    # fit L0 and L2 to the GOOD (reliable) bins only
    for L, col in [(0, "C0"), (2, "C2")]:
        a, c2v = fit_one(t[good], s[good], e[good], curves[(L, sidx)])
        ax.plot(theta_grid, a * curves[(L, sidx)](theta_grid), col, lw=2,
                label=f"L={L} fit (reliable bins): χ²/ν={c2v:.2f}")
    ax.errorbar(t[~bad], s[~bad], yerr=e[~bad], fmt="o", color="black",
                ms=7, capsize=3, label="reliable bins")
    ax.errorbar(t[bad], s[bad], yerr=e[bad], fmt="s", color="red",
                ms=9, capsize=3, label="13.75° bin (ε×4 swing)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
    ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")
    ax.set_title(f"{title}: L=0 vs L=2 once the 13.75° bin is excluded")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(HERE, "..", "plots", "eff_artifact_check.png")
fig.savefig(out, dpi=110)
print(f"saved: {out}")
