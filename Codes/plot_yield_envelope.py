#!/usr/bin/env python3
"""Per-state yield envelope across all systematic variants.

Reads ``results/<binning>[_<variant>]/yields.csv`` for every variant dir
matching ``<binning>*``, integrates each state's yield over all theta bins,
and renders:

    1. A compact summary plot ``plots/yield_envelope_<binning>.png``:
       nominal yield (square marker) per state, with the min..max envelope
       drawn as a gray bar and individual variants overlaid as small markers
       coloured by family (BW shape / PS / bg / L / cuts / contaminant).
    2. A CSV ``results/yield_envelope_<binning>.csv`` with one row per state:
       state, nominal, env_lo, env_hi, frac_lo, frac_hi, n_variants.
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
PLOTS = REPO / "plots"

# Family classification by variant name prefix/contents.  Order matters --
# the first matching rule wins.  "nominal" is a special case.
FAMILY_RULES = [
    ("nominal",                          re.compile(r"^$")),
    ("BW Γ scale",                  re.compile(r"\bbw_g\d")),
    ("phase space",                      re.compile(r"\bps[12]n_x")),
    ("linear bg scale",                  re.compile(r"\bbg_x")),
    ("drop 3.661 NEW",                   re.compile(r"\bno_3661\b")),
    ("L override",                       re.compile(r"\bL_s\d_L\d\b")),
]
FAMILY_COLORS = {
    "nominal":           "k",
    "BW Γ scale":   "#1f77b4",
    "phase space":       "#2ca02c",
    "linear bg scale":   "#9467bd",
    "drop 3.661 NEW":    "#ff7f0e",
    "L override":        "#8c564b",
    "strict cuts":       "#17becf",
    "contaminant":       "#d62728",
    "strict + contam":   "#e377c2",
}

ALL_STATES = [
    ("gs_3half+",        "gs 3/2$^+$"),
    ("ex2.763_1half-",   "2.763 1/2$^-$"),
    ("ex2.980_3half+",   "2.980 3/2$^+$"),
    ("ex3.661_NEW",      "3.661 NEW"),
    ("ex4.231_3half+",   "4.231 3/2$^+$"),
    ("ex4.841_1half+",   "4.841 1/2$^+$"),
    ("ex5.91_3half+",    "5.91 3/2$^+$"),
    ("ex6.30_5half+",    "6.30 5/2$^+$"),
]


def classify(variant_tag: str) -> str:
    """Map a variant tag ('' for nominal, otherwise the suffix) to a family."""
    # cuts / contam axis comes first; they dominate the prefix.
    if variant_tag.startswith("strict_contam"):  return "strict + contam"
    if variant_tag.startswith("contam"):         return "contaminant"
    if variant_tag.startswith("strict"):         return "strict cuts"
    for fam, rule in FAMILY_RULES:
        if rule.search(variant_tag):
            return fam
    return "other"


def state_totals(yields_csv: Path) -> dict[str, float]:
    """Sum integrated counts per state across all theta bins."""
    out = defaultdict(float)
    for row in csv.DictReader(yields_csv.open()):
        s = row["state"]
        if s.startswith("ex_"):  # bound 0.218 / 0.335 satellites
            continue
        out[s] += float(row["counts"])
    return dict(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", default="fine")
    ap.add_argument("--exclude-contam-state", action="store_true",
                    help="omit the 'contaminant' pseudo-state from the plot "
                         "(it is zero in non-contam variants)")
    # Default exclusions:
    #   * no_3661   -- this is a state-existence test, not a yield systematic.
    #                  Including it forces 3.661's envelope to 0 and inflates
    #                  4.231 by >100%.
    #   * contam_plus_relax -- known-unstable degeneracy between two knobs
    #                  trying to cover the same low-Ex excess; produces
    #                  unphysical jumps in 2.980 yield.
    ap.add_argument("--exclude", default="no_3661,contam_plus_relax",
                    help="comma list of variant-tag substrings to exclude "
                         "from the envelope (default keeps the envelope "
                         "clean against pathological variants)")
    args = ap.parse_args()
    exclude_subs = [s for s in args.exclude.split(",") if s]

    base = args.binning
    # Discover all variant dirs.  nominal sits at results/<base>/.
    # Variants sit at results/<base>_<tag>/ with tag the variant name.
    variants = {}
    for d in sorted(RESULTS.glob(f"{base}*")):
        if not d.is_dir():
            continue
        yc = d / "yields.csv"
        if not yc.exists():
            continue
        tag = "" if d.name == base else d.name[len(base) + 1:]
        if any(sub and sub in tag for sub in exclude_subs):
            print(f"  skip {d.name}  (matches exclusion)")
            continue
        variants[tag] = state_totals(yc)
    if not variants:
        raise SystemExit(f"no variant dirs found under results/{base}*")
    if "" not in variants:
        raise SystemExit(f"missing nominal: results/{base}/yields.csv")

    nominal = variants[""]
    states = list(ALL_STATES)
    if not args.exclude_contam_state:
        states = states + [("contaminant", "contaminant (0.45 MeV)")]

    # Envelope per state.
    print(f"{'state':<22} {'nominal':>9} {'env_lo':>9} {'env_hi':>9} "
          f"{'-%':>6} {'+%':>6}  n_variants={len(variants)}")
    summary_rows = []
    for key, _ in states:
        vals = [v.get(key, 0.0) for v in variants.values()]
        nom = nominal.get(key, 0.0)
        lo, hi = min(vals), max(vals)
        f_lo = (nom - lo) / nom * 100 if nom else 0
        f_hi = (hi - nom) / nom * 100 if nom else 0
        print(f"{key:<22} {nom:>9.1f} {lo:>9.1f} {hi:>9.1f} "
              f"{-f_lo:>+5.1f}% {+f_hi:>+5.1f}%")
        summary_rows.append({
            "state":       key,
            "nominal":     f"{nom:.2f}",
            "env_lo":      f"{lo:.2f}",
            "env_hi":      f"{hi:.2f}",
            "frac_lo_pct": f"{-f_lo:.2f}",
            "frac_hi_pct": f"{+f_hi:.2f}",
            "n_variants":  len(variants),
        })

    out_csv = RESULTS / f"yield_envelope_{base}.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader(); w.writerows(summary_rows)
    print(f"\nwrote {out_csv}")

    # ---------- plot ----------
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(states))
    label_seen = set()
    for i, (key, _) in enumerate(states):
        # envelope as a vertical bar
        vals = [v.get(key, 0.0) for v in variants.values()]
        lo, hi = min(vals), max(vals)
        ax.add_patch(plt.Rectangle((i - 0.30, lo), 0.60, hi - lo,
                                   facecolor="0.85", edgecolor="0.6",
                                   linewidth=0.5, zorder=1))
        # individual variants as small markers, jittered horizontally
        rng = np.random.default_rng(seed=i)
        for tag, vmap in variants.items():
            y = vmap.get(key, 0.0)
            fam = classify(tag)
            jx = i + rng.uniform(-0.20, 0.20)
            label = fam if fam not in label_seen else None
            ax.plot(jx, y, "o", ms=3, color=FAMILY_COLORS.get(fam, "0.5"),
                    alpha=0.7, label=label, zorder=3)
            if label: label_seen.add(label)
        # nominal as a black square on top
        ax.plot(i, nominal.get(key, 0.0), "s", ms=8,
                color="k", markerfacecolor="white",
                markeredgewidth=1.5, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in states], rotation=20, ha="right")
    ax.set_ylabel("Integrated counts (10$-$40 deg)")
    ax.set_title(f"Per-state yield envelope, {base} binning "
                 f"({len(variants)} systematic variants)")
    ax.grid(True, axis="y", alpha=0.3)
    handles, labels = ax.get_legend_handles_labels()
    # add nominal marker manually
    handles = [plt.Line2D([0], [0], marker="s", color="k", lw=0,
                          markerfacecolor="white", markeredgewidth=1.5,
                          ms=8, label="nominal"),
               plt.Rectangle((0, 0), 1, 1, facecolor="0.85", edgecolor="0.6",
                             label="envelope (min--max)")] + handles
    labels  = ["nominal", "envelope (min--max)"] + labels
    ax.legend(handles, labels, loc="upper left", fontsize=8, ncols=2)

    out_png = PLOTS / f"yield_envelope_{base}.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
