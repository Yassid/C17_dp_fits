"""
Spectral-fit quality comparison across systematic variants.

For each variant in run_systematics.VARIANTS, reads:
  - results/<binning>[_<variant>]/fit_quality.csv  (inclusive + per-bin chi2/NDF)
  - results/<binning>[_<variant>]/dsdo.csv        (extracted dsdo per state)

Produces:
  - results/fit_quality_summary_<binning>.csv    (one row per (variant, scope))
  - results/fit_quality_summary_<binning>.md     (compact table)
  - plots/syst_compare/fit_quality_<binning>.png (inclusive bar + per-bin lines)
  - plots/syst_compare/dsdo_spread_<binning>.png (per-state dsdo spread)
  - plots/syst_compare/dsdo_<state>_<binning>.png (per-state overlay)

Run from Codes/:
    python3 compare_variants_fitquality.py
    python3 compare_variants_fitquality.py --binning coarse2
"""
import argparse
import csv
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_systematics import VARIANTS, STATE_KEYS

# For a state explicitly dropped in a variant, its dsdo collapses to zero
# (efficiency*norm * tiny yield).  Skip those (variant, state) pairs in the
# spread computation -- they would otherwise dominate the band with what is
# really a "model alternative", not a "extracted-value uncertainty".
def _dropped_for_state(env_vars, state_idx):
    sd = env_vars.get("STATE_DROP", "")
    if not sd: return False
    parts = sd.split(",")
    if state_idx >= len(parts): return False
    try: return int(parts[state_idx]) != 0
    except ValueError: return False

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
RESULTS = REPO / "results"
PLOTS   = REPO / "plots" / "syst_compare"
PLOTS.mkdir(parents=True, exist_ok=True)

STATE_TITLES = {
    "gs_3half+":       "g.s. 3/2+",
    "ex_1half+_217":   "0.217 1/2+",
    "ex_5half+_335":   "0.335 5/2+",
    "ex2.763_1half-":  "2.763 1/2-",
    "ex2.980_3half+":  "2.980 3/2+",
    "ex3.661_NEW":     "3.661 NEW",
    "ex4.231_3half+":  "4.231 3/2+",
    "ex4.841_1half+":  "4.841 1/2+",
    "ex5.91_3half+":   "5.91 3/2+",
    "ex6.30_5half+":   "6.30 5/2+",
}


def variant_dir(binning, vname):
    sub = binning if vname == "nominal" else f"{binning}_{vname}"
    return RESULTS / sub


def load_fit_quality(binning, vname):
    p = variant_dir(binning, vname) / "fit_quality.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def load_dsdo(binning, vname):
    p = variant_dir(binning, vname) / "dsdo.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", default="fine")
    args = ap.parse_args()
    bn = args.binning

    # ---- collect ----
    names = [n for n, _ in VARIANTS]
    fq = {n: load_fit_quality(bn, n) for n in names}
    dd = {n: load_dsdo(bn, n)        for n in names}
    missing = [n for n in names if fq[n] is None]
    if missing:
        print(f"WARNING: missing fit_quality.csv for: {missing}")
    names = [n for n in names if fq[n] is not None]

    # ---- summary CSV ----
    rows = []
    for n in names:
        df = fq[n]
        incl = df[df["scope"] == "inclusive"].iloc[0]
        bins = df[df["scope"] == "bin"]
        rows.append({
            "variant": n,
            "inclusive_chi2": float(incl["chi2"]),
            "inclusive_ndf":  int(incl["ndf"]),
            "inclusive_chi2_per_ndf": float(incl["chi2_per_ndf"]),
            "perbin_median_chi2_per_ndf": float(bins["chi2_per_ndf"].median()),
            "perbin_max_chi2_per_ndf":    float(bins["chi2_per_ndf"].max()),
            "perbin_min_chi2_per_ndf":    float(bins["chi2_per_ndf"].min()),
            "total_chi2":   float(bins["chi2"].sum() + incl["chi2"]),
            "total_ndf":    int(bins["ndf"].sum() + incl["ndf"]),
        })
    sumdf = pd.DataFrame(rows)
    sumdf["total_chi2_per_ndf"] = sumdf["total_chi2"] / sumdf["total_ndf"]

    out_csv = RESULTS / f"fit_quality_summary_{bn}.csv"
    sumdf.to_csv(out_csv, index=False, float_format="%.4f")
    print(f"wrote {out_csv}")

    # ---- markdown table (vs nominal) ----
    nom = sumdf[sumdf["variant"] == "nominal"].iloc[0]
    out_md = RESULTS / f"fit_quality_summary_{bn}.md"
    with out_md.open("w") as fh:
        fh.write(f"# Spectral-fit quality across variants ({bn})\n\n")
        fh.write("Columns: inclusive chi2/NDF, per-bin chi2/NDF (median, max), "
                 "and total chi2/NDF summed over inclusive + all theta_CM bins. "
                 "Delta(total chi2/NDF) is relative to nominal.\n\n")
        fh.write("| variant | inclusive chi2/NDF | bin chi2/NDF (med, max) | "
                 "total chi2/NDF | Delta vs nominal |\n")
        fh.write("|---|---|---|---|---|\n")
        for _, r in sumdf.iterrows():
            dlt = r["total_chi2_per_ndf"] - nom["total_chi2_per_ndf"]
            mark = ""
            if abs(dlt) > 0.10: mark = " **"
            fh.write(f"| {r['variant']} | "
                     f"{r['inclusive_chi2']:.1f} / {r['inclusive_ndf']} = {r['inclusive_chi2_per_ndf']:.2f} | "
                     f"{r['perbin_median_chi2_per_ndf']:.2f}, {r['perbin_max_chi2_per_ndf']:.2f} | "
                     f"{r['total_chi2_per_ndf']:.3f} | "
                     f"{dlt:+.3f}{mark} |\n")
    print(f"wrote {out_md}")

    # ---- plot: inclusive chi2/NDF bar + per-bin lines ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9),
                                   gridspec_kw=dict(height_ratios=[1, 2]))
    # ax1: bar of inclusive chi2/NDF per variant (color = vs nominal)
    nom_inc = nom["inclusive_chi2_per_ndf"]
    colors = ["C7" if v == "nominal"
              else "C0" if (sumdf.loc[sumdf["variant"]==v, "inclusive_chi2_per_ndf"].values[0] < nom_inc + 0.10)
              else "C3"
              for v in names]
    ax1.bar(names, sumdf["inclusive_chi2_per_ndf"], color=colors, edgecolor="k", lw=0.4)
    ax1.axhline(nom_inc, color="k", ls="--", alpha=0.5, label=f"nominal {nom_inc:.2f}")
    ax1.set_ylabel("inclusive chi2/NDF")
    ax1.set_title(f"Inclusive Ex-spectrum fit quality across variants ({bn})")
    ax1.tick_params(axis="x", labelrotation=70, labelsize=8)
    ax1.legend(fontsize=9)
    ax1.grid(axis="y", alpha=0.3)

    # ax2: per-bin chi2/NDF vs theta_CM, one curve per variant
    cmap = plt.get_cmap("tab20")
    for i, n in enumerate(names):
        df = fq[n]
        bdf = df[df["scope"] == "bin"].sort_values("bin_center")
        style = dict(marker="o", ms=4, lw=1.0, alpha=0.7)
        if n == "nominal":
            style.update(color="k", lw=2.2, alpha=1.0, ms=6, zorder=5,
                         label="nominal")
        else:
            style.update(color=cmap(i % 20), label=n)
        ax2.plot(bdf["bin_center"], bdf["chi2_per_ndf"], **style)
    ax2.axhline(1.0, color="grey", ls=":", alpha=0.5)
    ax2.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
    ax2.set_ylabel("per-bin chi2/NDF")
    ax2.set_title("Per-bin spectral fit quality across variants")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=6.5, ncol=4, loc="upper left")
    fig.tight_layout()
    out_png = PLOTS / f"fit_quality_{bn}.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"saved: {out_png}")

    # ---- per-state dsdo overlay across variants ----
    states_to_plot = [s for s in STATE_KEYS if dd["nominal"] is not None
                      and (dd["nominal"]["state"] == s).any()]
    n_states = len(states_to_plot)
    ncol = 4
    nrow = (n_states + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5*ncol, 3.6*nrow), sharex=True)
    axes = np.atleast_2d(axes)
    band_rows = []   # for systematic-band CSV
    # Map variant name -> env-vars dict (for dropped-state filtering)
    variant_env = {n: ev for n, ev in VARIANTS}
    for i, state in enumerate(states_to_plot):
        ax = axes[i // ncol, i % ncol]
        nom_df = dd["nominal"][dd["nominal"]["state"] == state].sort_values("bin_center")
        theta_nom = nom_df["bin_center"].to_numpy()
        s_nom     = nom_df["dsdo_mbsr"].to_numpy()
        e_nom     = nom_df["dsdo_err_mbsr"].to_numpy()
        # state_idx in the 7-unbound list, for dropped-state filtering
        state_idx = STATE_KEYS.index(state) if state in STATE_KEYS else -1
        # collect every variant at the same theta points
        variant_grid = []
        for n in names:
            df = dd[n]
            if df is None: continue
            # skip variant where this very state was dropped
            if state_idx >= 0 and _dropped_for_state(variant_env.get(n, {}), state_idx):
                continue
            sub = df[df["state"] == state].sort_values("bin_center")
            # align on theta
            row = []
            for t in theta_nom:
                m = (sub["bin_center"] == t)
                row.append(float(sub.loc[m, "dsdo_mbsr"].iloc[0]) if m.any() else np.nan)
            variant_grid.append((n, np.array(row)))
        # plot variants light
        for n, row in variant_grid:
            if n == "nominal": continue
            ax.plot(theta_nom, row, color="C0", lw=0.6, alpha=0.35)
        # systematic band: min/max across variants
        stack = np.vstack([row for _, row in variant_grid])
        lo = np.nanmin(stack, axis=0)
        hi = np.nanmax(stack, axis=0)
        ax.fill_between(theta_nom, lo, hi, color="C0", alpha=0.18,
                        label="syst spread")
        # nominal with stat error
        ax.errorbar(theta_nom, s_nom, yerr=e_nom, fmt="o", color="k",
                    ms=5, capsize=2, label="nominal (stat)")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 50)
        ax.set_xlim(8, 42)
        ax.set_title(STATE_TITLES.get(state, state), fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
        if i // ncol == nrow - 1:
            ax.set_xlabel(r"$\theta_\mathrm{CM}$ (deg)")
        if i % ncol == 0:
            ax.set_ylabel(r"$d\sigma/d\Omega$ (mb/sr)")
        if i == 0:
            ax.legend(fontsize=7.5, loc="lower left")
        # band CSV row
        for k, t in enumerate(theta_nom):
            stat_err = e_nom[k]
            syst_half = 0.5 * (hi[k] - lo[k]) if np.isfinite(hi[k] - lo[k]) else float("nan")
            band_rows.append(dict(
                state=state,
                bin_center=t,
                nominal=s_nom[k],
                stat_err=stat_err,
                syst_half_range=syst_half,
                syst_lo=lo[k], syst_hi=hi[k],
                syst_over_stat=(syst_half/stat_err if stat_err > 0 else float("nan")),
            ))
    # hide empties
    for j in range(n_states, nrow*ncol):
        axes[j // ncol, j % ncol].axis("off")
    fig.suptitle(f"Per-state dsdo: nominal (stat) vs systematic spread across {len(names)-1} variants ({bn})",
                 fontsize=12)
    fig.tight_layout()
    out_png = PLOTS / f"dsdo_spread_{bn}.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"saved: {out_png}")

    band = pd.DataFrame(band_rows)
    out_csv = RESULTS / f"systematics_dsdo_band_{bn}.csv"
    band.to_csv(out_csv, index=False, float_format="%.6g")
    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
