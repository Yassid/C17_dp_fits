"""
Overlay FRESCO ADWA shapes (weakly-bound approximation, be=0.5 MeV) for the 7
17C unbound resonances against the per-state angular distributions extracted
from C16_pd_AngBins (dsdo.csv).  For each state, fit a single normalization
factor (effective C²S) to the data.

FRESCO states 1..7 map to:
  1: 2.763 1/2-   (2p1/2, L=1)
  2: 2.980 3/2+   (2d3/2, L=2)
  3: 3.661 NEW    (2d3/2, L=2)
  4: 4.231 3/2+   (2d3/2, L=2)
  5: 4.841 1/2+   (3s1/2, L=0)
  6: 5.91  3/2+   (2d3/2, L=2)
  7: 6.30  5/2+   (2d5/2, L=2)
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

state_names = [
    "ex2.763_1half-",
    "ex2.980_3half+",
    "ex3.661_NEW",
    "ex4.231_3half+",
    "ex4.841_1half+",
    "ex5.91_3half+",
    "ex6.30_5half+",
]
state_titles = [
    r"2.763 1/2$^-$  (2p1/2, L=1)",
    r"2.980 3/2$^+$  (2d3/2, L=2)",
    r"3.661 NEW  (2d3/2, L=2)",
    r"4.231 3/2$^+$  (2d3/2, L=2)",
    r"4.841 1/2$^+$  (3s1/2, L=0)",
    r"5.91 3/2$^+$  (2d3/2, L=2)",
    r"6.30 5/2$^+$  (2d5/2, L=2)",
]

df = pd.read_csv(os.path.join(RESULTS_DIR, "dsdo.csv"))

# Load FRESCO curves
curves = []
for i in range(1, 8):
    th, sig = np.loadtxt(os.path.join(FRESCO, f"dp17C_unbound_adwa_state{i}.dat")).T
    curves.append(interp1d(th, sig, bounds_error=False, fill_value=np.nan))

# Fit one normalization per state, restricting to forward-half bins where the
# AT-TPC kinematic edge does not bias the data toward the bg floor.
# (32.5+ bins for unbound have amplitudes near floor.)
theta_grid = np.linspace(1, 50, 600)

def fit_norm(theta, sig, err, curve, kmin=None, kmax=None):
    # mask out NaN and the kinematic-edge floor
    floor = 0.03   # mb/sr -- below this is at the bin-amplitude floor
    mask = (sig > floor) & (err > 0)
    if kmin is not None: mask &= (theta >= kmin)
    if kmax is not None: mask &= (theta <= kmax)
    if mask.sum() < 2:
        return 0.0, 0.0, float("nan"), 0
    s = sig[mask]; e = err[mask]; t = theta[mask]
    m = curve(t)
    W = 1.0 / e ** 2
    num = np.nansum(s * m * W)
    den = np.nansum(m * m * W)
    a = num / den if den > 0 else 0.0
    chi2 = float(np.nansum(((s - a * m) / e) ** 2))
    ndf = int(mask.sum() - 1)
    da = float(np.sqrt(1.0 / den)) if den > 0 else 0.0
    return a, da, chi2 / max(ndf, 1), ndf


fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=False)

print("FRESCO ADWA weakly-bound (be=0.5 MeV) fits to per-state AD")
print(f"{'state':<28s} {'kmax':>6s} {'C²S_eff':>10s} {'±':>10s} {'χ²/ν':>8s} {'ν':>4s}")
for i, (sname, stitle) in enumerate(zip(state_names, state_titles)):
    ax = axes[i // 4, i % 4]
    sub = df[df["state"] == sname].sort_values("bin_center")
    theta = sub["bin_center"].to_numpy()
    sig   = sub["dsdo_mbsr"].to_numpy()
    err   = sub["dsdo_err_mbsr"].to_numpy()

    a, da, c2v, ndf = fit_norm(theta, sig, err, curves[i], kmax=30.0)

    ax.errorbar(theta, sig, yerr=err, fmt="o", color="black", ms=6, capsize=2,
                label="data (12-bin fit)")
    if np.isfinite(a) and a > 0:
        ax.plot(theta_grid, a * curves[i](theta_grid), "r-", lw=2,
                label=(f"FRESCO  C²S_eff = {a:.3f} ± {da:.3f}\n"
                       f"χ²/ν = {c2v:.2f}"))
    ax.plot(theta_grid, curves[i](theta_grid), "r--", lw=1, alpha=0.4,
            label="FRESCO (C²S=1)")
    ax.set_yscale("log")
    ax.set_xlim(8, 42)
    ax.set_ylim(1e-3, 1e3)
    ax.set_title(stitle, fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    if i // 4 == 1: ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
    if i % 4 == 0:  ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")
    ax.legend(fontsize=7.5, loc="upper right")

    print(f"{sname:<28s} {30.0:>6.1f} {a:>10.4f} {da:>10.4f} {c2v:>8.2f} {ndf:>4d}")

# Last panel: pure FRESCO shapes overlaid, normalized to peak=1
ax = axes[1, 3]
for i, stitle in enumerate(state_titles):
    y = curves[i](theta_grid)
    y_norm = y / np.nanmax(y)
    ax.plot(theta_grid, y_norm, lw=1.7, label=stitle.split()[0] + " " + stitle.split()[1])
ax.set_yscale("log"); ax.set_xlim(0, 50); ax.set_ylim(5e-3, 1.3)
ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)"); ax.set_ylabel("dσ/dΩ (peak-normalised)")
ax.set_title("FRESCO shapes (peak=1)", fontsize=10)
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)

fig.suptitle("16C(d,p)17C @ 11.77 MeV/u — FRESCO ADWA (weakly-bound approx) vs extracted AD",
             fontsize=12)
fig.tight_layout()
out = os.path.join(PLOTS_DIR, "plots_fresco_unbound.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nsaved: {out}")
