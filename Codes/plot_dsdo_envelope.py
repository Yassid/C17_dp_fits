#!/usr/bin/env python3
"""Per-state dsdo angular distribution with systematic envelope band.

Reads ``results/<binning>[_<variant>]/dsdo.csv`` across all variant dirs
and renders a multi-panel plot (one panel per unbound state) showing:
  - nominal points with statistical error bars
  - shaded envelope band from min..max across all variants at each theta bin

Output: ``plots/dsdo_envelope_<binning>.png``
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
PLOTS = REPO / "plots"

STATES = [
    ("ex2.763_1half-",   "2.763 1/2$^-$"),
    ("ex2.980_3half+",   "2.980 3/2$^+$"),
    ("ex3.661_NEW",      "3.661 NEW"),
    ("ex4.231_3half+",   "4.231 3/2$^+$"),
    ("ex4.841_1half+",   "4.841 1/2$^+$"),
    ("ex5.91_3half+",    "5.91 3/2$^+$"),
    ("ex6.30_5half+",    "6.30 5/2$^+$"),
]


def load(dirpath: Path) -> dict[str, list[tuple[float, float, float]]]:
    """Per state -> list of (theta_center, dsdo, dsdo_err)."""
    out = defaultdict(list)
    fp = dirpath / "dsdo.csv"
    if not fp.exists():
        return {}
    for r in csv.DictReader(fp.open()):
        s = r["state"]
        try:
            th = float(r["bin_center"])
            ds = float(r["dsdo_mbsr"])
            de = float(r["dsdo_err_mbsr"])
        except (KeyError, ValueError):
            continue
        out[s].append((th, ds, de))
    # Sort by theta for safety.
    return {s: sorted(v) for s, v in out.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", default="fine")
    ap.add_argument("--exclude", default="no_3661,contam_plus_relax",
                    help="comma list of variant-tag substrings to exclude "
                         "(same default as plot_yield_envelope.py)")
    args = ap.parse_args()
    exclude_subs = [s for s in args.exclude.split(",") if s]

    base = args.binning
    variants = {}
    for d in sorted(RESULTS.glob(f"{base}*")):
        if not d.is_dir():
            continue
        rows = load(d)
        if not rows:
            continue
        tag = "" if d.name == base else d.name[len(base) + 1:]
        if any(sub and sub in tag for sub in exclude_subs):
            continue
        variants[tag] = rows
    if "" not in variants:
        raise SystemExit(f"missing nominal: results/{base}/dsdo.csv")

    nominal = variants[""]

    # For envelope we need a common theta grid.  Use nominal's grid and
    # interpolate / look up other variants at the same theta_center (the
    # binning is identical across variants, so theta centers match).
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), sharex=True)
    axes = axes.flatten()
    for ax, (key, label) in zip(axes, STATES):
        if key not in nominal:
            ax.set_visible(False); continue
        nom = nominal[key]
        th_nom = np.array([t for t, _, _ in nom])
        ds_nom = np.array([d for _, d, _ in nom])
        de_nom = np.array([e for _, _, e in nom])

        # Gather all variants at each nominal theta.
        env_lo = np.full_like(ds_nom, np.inf)
        env_hi = np.full_like(ds_nom, -np.inf)
        for tag, vmap in variants.items():
            v = vmap.get(key)
            if not v:
                continue
            v_th = np.array([t for t, _, _ in v])
            v_ds = np.array([d for _, d, _ in v])
            # match thetas exactly (variant dropouts produce gaps)
            for i, t in enumerate(th_nom):
                idx = np.where(np.isclose(v_th, t, atol=0.05))[0]
                if not len(idx):
                    continue
                y = v_ds[idx[0]]
                env_lo[i] = min(env_lo[i], y)
                env_hi[i] = max(env_hi[i], y)
        # Replace inf with nominal in case a theta had no variants
        bad = ~np.isfinite(env_lo)
        env_lo[bad] = ds_nom[bad]
        env_hi[bad] = ds_nom[bad]

        ax.fill_between(th_nom, env_lo, env_hi, color="0.80",
                        label="syst envelope", zorder=1)
        ax.errorbar(th_nom, ds_nom, yerr=de_nom, fmt="o", ms=4,
                    color="k", capsize=2, label="nominal $\\pm$ stat",
                    zorder=3)
        ax.set_yscale("log")
        ax.set_title(label, fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
        if ax is axes[0]:
            ax.legend(fontsize=8, loc="upper right")
        if ax in axes[4:]:
            ax.set_xlabel(r"$\theta_{\rm CM}$ (deg)")
        if ax in (axes[0], axes[4]):
            ax.set_ylabel(r"d$\sigma$/d$\Omega$ (mb/sr)")

    # Hide unused subplot
    if len(STATES) < len(axes):
        for ax in axes[len(STATES):]:
            ax.set_visible(False)

    fig.suptitle(f"Angular distributions with systematic envelope "
                 f"(n_variants = {len(variants)}, binning = {base})",
                 y=1.00)
    fig.tight_layout()
    out = PLOTS / f"dsdo_envelope_{base}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
