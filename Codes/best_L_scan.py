"""
For each of the 7 unbound 17C resonances, fit the extracted angular
distribution to FRESCO ADWA shapes for L = 0, 1, 2, 3 (weakly-bound approx).
Report the best L by χ²/ν and produce a side-by-side plot.

Inputs:
  - /home/yassid/C16dp_fits/Codes/dsdo.csv          (extracted data)
  - /home/yassid/fair_install/16C_dp/fresco/outputs/dp17C_L{L}_state{S}.dat
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
BINNING = os.environ.get("BINNING", "fine")
RESULTS_DIR = os.path.join(HERE, "..", "results", BINNING)
PLOTS_DIR = os.path.join(HERE, "..", "plots", BINNING)
os.makedirs(PLOTS_DIR, exist_ok=True)

state_keys = [
    "ex2.763_1half-",   # 0
    "ex2.980_3half+",
    "ex3.661_NEW",
    "ex4.231_3half+",
    "ex4.841_1half+",
    "ex5.91_3half+",
    "ex6.30_5half+",
]
state_titles = [
    r"2.763 1/2$^-$",
    r"2.980 3/2$^+$",
    r"3.661 NEW",
    r"4.231 3/2$^+$",
    r"4.841 1/2$^+$",
    r"5.91 3/2$^+$",
    r"6.30 5/2$^+$",
]
L_labels = ["L=0 (3s1/2)", "L=1 (2p1/2)", "L=2 (2d3/2)", "L=3 (1f5/2)"]
L_colors = ["C0", "C1", "C2", "C3"]
KMAX = 30.0      # fit only forward bins -- avoid kinematic-edge floor
TIE_THRESH = 0.3  # if the runner-up L is within this in reduced chi2, the
                  # L assignment is a tie (shapes statistically indistinct)

df = pd.read_csv(os.path.join(RESULTS_DIR, "dsdo.csv"))

# Load all (L, state) FRESCO curves
curves = {}
for L in range(4):
    for s in range(1, 8):
        th, sig = np.loadtxt(os.path.join(FRESCO, f"dp17C_L{L}_state{s}.dat")).T
        curves[(L, s)] = interp1d(th, sig, bounds_error=False, fill_value=np.nan)


def fit_one(theta, sig, err, curve):
    floor = 0.03
    mask = (sig > floor) & (err > 0) & (theta <= KMAX)
    if mask.sum() < 2:
        return 0.0, 0.0, float("nan"), 0
    s, e, t = sig[mask], err[mask], theta[mask]
    m = curve(t)
    W = 1.0 / e ** 2
    num = float(np.nansum(s * m * W))
    den = float(np.nansum(m * m * W))
    if den <= 0:
        return 0.0, 0.0, float("nan"), 0
    a = num / den
    chi2 = float(np.nansum(((s - a * m) / e) ** 2))
    ndf = int(mask.sum() - 1)
    da = float(np.sqrt(1.0 / den))
    return a, da, chi2 / max(ndf, 1), ndf


theta_grid = np.linspace(1, 50, 600)

print("=" * 80)
print(f"Per-state L scan  (data points kept: θ_CM ≤ {KMAX:.0f}°, σ > 0.03 mb/sr floor)")
print(f"{'state':<22s} | " + " | ".join(f"{L_labels[L]:>13s} χ²/ν" for L in range(4)) + " | best")
print("-" * 80)

# Side-by-side per-state plot
fig, axes = plt.subplots(2, 4, figsize=(20, 9))
results = []
for i, (key, title) in enumerate(zip(state_keys, state_titles)):
    ax = axes[i // 4, i % 4]
    sub = df[df["state"] == key].sort_values("bin_center")
    theta = sub["bin_center"].to_numpy()
    sig   = sub["dsdo_mbsr"].to_numpy()
    err   = sub["dsdo_err_mbsr"].to_numpy()
    ax.errorbar(theta, sig, yerr=err, fmt="o", color="black", ms=6, capsize=2,
                label="data")

    fits = []
    for L in range(4):
        a, da, c2v, ndf = fit_one(theta, sig, err, curves[(L, i + 1)])
        fits.append((L, a, da, c2v, ndf))
        if np.isfinite(c2v) and a > 0:
            ax.plot(theta_grid, a * curves[(L, i + 1)](theta_grid),
                    color=L_colors[L], lw=1.6,
                    label=f"{L_labels[L]}: C²S={a:.3g}, χ²/ν={c2v:.2f}")

    # Best L = smallest chi2 with positive amp.  Flag a tie when the
    # runner-up is within TIE_THRESH in reduced chi2: the L shapes are
    # then statistically indistinguishable for this state, so a single
    # "best L" would be misleading (see eff_artifact_check.py).
    valid = sorted([f for f in fits if np.isfinite(f[3]) and f[1] > 0],
                   key=lambda x: x[3])
    if valid:
        bestL = valid[0][0]
        tie = len(valid) > 1 and (valid[1][3] - valid[0][3]) < TIE_THRESH
        tieL = valid[1][0] if tie else -1
    else:
        bestL, tie, tieL = -1, False, -1
    results.append((title, fits, bestL, tie, tieL))

    if bestL < 0:
        Llabel = "none"
    elif tie:
        Llabel = f"L={bestL}/{tieL} (tie)"
    else:
        Llabel = f"L={bestL}"

    ax.set_yscale("log"); ax.set_xlim(8, 42); ax.set_ylim(1e-3, 5e2)
    ax.set_title(f"{title}   {Llabel}", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.0, loc="upper right")
    if i // 4 == 1: ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
    if i % 4 == 0:  ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")

    row = f"{key:<22s} | "
    for L, a, da, c2v, ndf in fits:
        cell = f"{c2v:6.2f}" if np.isfinite(c2v) else "  nan"
        marker = " *" if L == bestL else (" ~" if L == tieL else "  ")
        row += f"{a:7.3g} {cell}{marker}| "
    row += f" {Llabel}"
    print(row)

axes[1, 3].axis("off")
axes[1, 3].text(0.0, 0.95, "Best-L summary", fontsize=13, weight="bold", transform=axes[1, 3].transAxes)
txt = "state           L assignment      χ²/ν\n"
txt += "-" * 56 + "\n"
for (title, fits, bestL, tie, tieL) in results:
    if bestL < 0:
        txt += f"{title:<14s}  none\n"
        continue
    b = fits[bestL]
    if tie:
        t = fits[tieL]
        txt += f"{title:<14s}  L={bestL}/{tieL} (tie)   {b[3]:.2f} / {t[3]:.2f}\n"
    else:
        txt += f"{title:<14s}  L={bestL}            {b[3]:.2f}\n"
axes[1, 3].text(0.0, 0.05, txt, fontsize=9.5, family="monospace",
                transform=axes[1, 3].transAxes, va="bottom")

fig.suptitle("16C(d,p)17C — L-scan per unbound state (FRESCO ADWA weakly-bound, be=0.5 MeV)",
             fontsize=12)
fig.tight_layout()
out = os.path.join(PLOTS_DIR, "plots_L_scan.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nsaved: {out}")

# Machine-readable summary CSV, ingested by the systematics driver.
csv_path = os.path.join(RESULTS_DIR, "L_scan.csv")
with open(csv_path, "w") as fh:
    fh.write("state,bestL,tie,tieL,c2v_L0,c2v_L1,c2v_L2,c2v_L3,amp_L0,amp_L1,amp_L2,amp_L3\n")
    for key, (title, fits, bestL, tie, tieL) in zip(state_keys, results):
        c2vs = ["nan", "nan", "nan", "nan"]
        amps = ["0",   "0",   "0",   "0"]
        for L, a, da, c2v, ndf in fits:
            c2vs[L] = f"{c2v:.4f}" if np.isfinite(c2v) else "nan"
            amps[L] = f"{a:.6g}"
        fh.write(f"{key},{bestL},{int(tie)},{tieL}," + ",".join(c2vs) + "," + ",".join(amps) + "\n")
print(f"saved: {csv_path}")
