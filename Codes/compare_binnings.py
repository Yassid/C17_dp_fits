"""
Side-by-side comparison of the dsdo angular distributions and FRESCO
single-norm fits across the three theta_CM binning schemes:
  fine     12 x 2.5 deg
  coarse2   6 x 5.0 deg
  coarse4   3 x 10  deg
Reads ../results/<scheme>/dsdo.csv produced by C16_pd_AngBins.C(scheme).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
FRESCO = os.path.join(HERE, "..", "fresco", "outputs")
SCHEMES = ["fine", "coarse2", "coarse4"]
SCHEME_LABELS = {
    "fine":    r"fine 12 $\times$ 2.5$^\circ$",
    "coarse2": r"coarse2 6 $\times$ 5$^\circ$",
    "coarse4": r"coarse4 3 $\times$ 10$^\circ$",
}
SCHEME_COLOR  = {"fine": "C0", "coarse2": "C1", "coarse4": "C3"}
SCHEME_MARKER = {"fine": "o",  "coarse2": "s",  "coarse4": "^"}

state_keys = [
    "gs_3half+", "ex_1half+_217", "ex_5half+_335",
    "ex2.763_1half-", "ex2.980_3half+", "ex3.661_NEW",
    "ex4.231_3half+", "ex4.841_1half+", "ex5.91_3half+", "ex6.30_5half+",
]
state_titles = {
    "gs_3half+":       r"g.s. 3/2$^+$",
    "ex_1half+_217":   r"0.217 1/2$^+$",
    "ex_5half+_335":   r"0.335 5/2$^+$",
    "ex2.763_1half-":  r"2.763 1/2$^-$  (L=1)",
    "ex2.980_3half+":  r"2.980 3/2$^+$",
    "ex3.661_NEW":     r"3.661 NEW",
    "ex4.231_3half+":  r"4.231 3/2$^+$  (L=2)",
    "ex4.841_1half+":  r"4.841 1/2$^+$  (L=0)",
    "ex5.91_3half+":   r"5.91 3/2$^+$",
    "ex6.30_5half+":   r"6.30 5/2$^+$",
}

# FRESCO unbound state files mapped to the 7 BW states (same order as in
# fresco_overlay_unbound.py).  Used to overlay the absolute-norm fit shape.
unbound_keys = [
    "ex2.763_1half-", "ex2.980_3half+", "ex3.661_NEW", "ex4.231_3half+",
    "ex4.841_1half+", "ex5.91_3half+", "ex6.30_5half+",
]
fresco_curves = {}
for i, key in enumerate(unbound_keys, start=1):
    fname = os.path.join(FRESCO, f"dp17C_unbound_adwa_state{i}.dat")
    if os.path.exists(fname):
        th, sig = np.loadtxt(fname).T
        fresco_curves[key] = interp1d(th, sig, bounds_error=False, fill_value=np.nan)

# Load per-scheme dsdo tables.
dfs = {s: pd.read_csv(os.path.join(HERE, "..", "results", s, "dsdo.csv"))
       for s in SCHEMES}


def fit_norm(theta, sig, err, curve, kmax=30.0, floor=0.03):
    mask = (sig > floor) & (err > 0) & (theta <= kmax)
    if mask.sum() < 2:
        return float("nan"), float("nan"), float("nan"), 0
    s, e, t = sig[mask], err[mask], theta[mask]
    m = curve(t)
    valid = ~np.isnan(m)
    if valid.sum() < 2:
        return float("nan"), float("nan"), float("nan"), 0
    s, e, t, m = s[valid], e[valid], t[valid], m[valid]
    W = 1.0 / e ** 2
    num = float(np.sum(s * m * W))
    den = float(np.sum(m * m * W))
    if den <= 0:
        return float("nan"), float("nan"), float("nan"), 0
    a = num / den
    chi2 = float(np.sum(((s - a * m) / e) ** 2))
    ndf = int(len(s) - 1)
    da = float(np.sqrt(1.0 / den))
    return a, da, chi2 / max(ndf, 1), ndf


theta_grid = np.linspace(1, 50, 600)

# 2x5 grid: one panel per state.
fig, axes = plt.subplots(2, 5, figsize=(22, 9), sharex=True)
table_rows = []

for i, state in enumerate(state_keys):
    ax = axes[i // 5, i % 5]
    for sc in SCHEMES:
        sub = dfs[sc][dfs[sc]["state"] == state].sort_values("bin_center")
        if sub.empty:
            continue
        ax.errorbar(sub["bin_center"], sub["dsdo_mbsr"], yerr=sub["dsdo_err_mbsr"],
                    fmt=SCHEME_MARKER[sc], color=SCHEME_COLOR[sc], ms=7, capsize=2.5,
                    label=SCHEME_LABELS[sc], alpha=0.85)
    # FRESCO overlay for unbound states: refit single norm per scheme.
    if state in fresco_curves:
        curve = fresco_curves[state]
        ax.plot(theta_grid, curve(theta_grid), color="k", lw=0.8, ls=":",
                alpha=0.4, label="FRESCO (C²S=1)")
        for sc in SCHEMES:
            sub = dfs[sc][dfs[sc]["state"] == state].sort_values("bin_center")
            a, da, c2v, ndf = fit_norm(sub["bin_center"].to_numpy(),
                                       sub["dsdo_mbsr"].to_numpy(),
                                       sub["dsdo_err_mbsr"].to_numpy(),
                                       curve)
            table_rows.append((state, sc, a, da, c2v, ndf))
            if np.isfinite(a) and a > 0:
                ax.plot(theta_grid, a * curve(theta_grid),
                        color=SCHEME_COLOR[sc], lw=1.2, alpha=0.6)
    ax.set_yscale("log")
    ax.set_xlim(8, 42)
    ax.set_ylim(1e-3, 50)
    ax.set_title(state_titles.get(state, state), fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    if i % 5 == 0:
        ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")
    if i // 5 == 1:
        ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
    if i == 0:
        ax.legend(fontsize=7.5, loc="lower left")

fig.suptitle(
    "16C(d,p)17C — dsigma/dOmega per state across three theta_CM binning schemes\n"
    "(thin coloured lines = FRESCO ADWA single-norm fit per scheme; dotted = FRESCO C^2S=1)",
    fontsize=12,
)
fig.tight_layout()
out_dir = os.path.join(HERE, "..", "plots", "comparison")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "compare_binnings.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"saved: {out}")

# --- summary table -------------------------------------------------
print("\nC²S_eff (FRESCO single-norm fit, θ_CM ≤ 30°) per state per scheme:")
print(f"{'state':<22s} {'scheme':<9s} {'C²S_eff':>10s} {'±':>10s} {'χ²/ν':>8s} {'ν':>4s}")
for state, sc, a, da, c2v, ndf in table_rows:
    a_s   = f"{a:>10.4f}"  if np.isfinite(a)  else f"{'-':>10s}"
    da_s  = f"{da:>10.4f}" if np.isfinite(da) else f"{'-':>10s}"
    c2_s  = f"{c2v:>8.2f}" if np.isfinite(c2v) else f"{'-':>8s}"
    print(f"{state:<22s} {sc:<9s} {a_s} {da_s} {c2_s} {ndf:>4d}")

# Save a tidy CSV as well.
tab = pd.DataFrame(table_rows, columns=["state", "scheme", "C2S_eff", "C2S_err", "chi2_red", "ndf"])
csv_out = os.path.join(out_dir, "fresco_C2Seff_compare.csv")
tab.to_csv(csv_out, index=False, float_format="%.6g")
print(f"\nsaved: {csv_out}")
