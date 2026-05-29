#!/usr/bin/env python3
"""Build a slides-ready Markdown summary of the spectral-fit systematics.

Reads:
  - results/yield_envelope_<binning>.csv (from plot_yield_envelope.py)
  - results/<binning>[_<variant>]/fit_quality.csv (per-variant inclusive chi2)
Writes:
  - results/systematics_summary_envelope_<binning>.md
"""
import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

STATE_LABELS = {
    "gs_3half+":        "g.s. 3/2+",
    "ex2.763_1half-":   "2.763 1/2-",
    "ex2.980_3half+":   "2.980 3/2+",
    "ex3.661_NEW":      "3.661 NEW",
    "ex4.231_3half+":   "4.231 3/2+",
    "ex4.841_1half+":   "4.841 1/2+",
    "ex5.91_3half+":    "5.91 3/2+",
    "ex6.30_5half+":    "6.30 5/2+",
    "contaminant":      "contaminant (0.45 MeV)",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", default="fine")
    args = ap.parse_args()
    base = args.binning

    env_csv = RESULTS / f"yield_envelope_{base}.csv"
    if not env_csv.exists():
        raise SystemExit(f"missing {env_csv} -- run plot_yield_envelope.py first")
    env = list(csv.DictReader(env_csv.open()))

    # Headline chi2 grid: nominal cuts vs strict, no-contam vs +contam.
    cells = {}
    for cut_tag in ("", "strict"):
        for cm_tag in ("", "contam"):
            parts = [t for t in [cut_tag, cm_tag] if t]
            name = "_".join(parts) if parts else ""
            d = RESULTS / (base if not name else f"{base}_{name}")
            fq = d / "fit_quality.csv"
            if not fq.exists():
                cells[(cut_tag, cm_tag)] = None
                continue
            inc = [r for r in csv.DictReader(fq.open()) if r["scope"] == "inclusive"][0]
            cells[(cut_tag, cm_tag)] = float(inc["chi2_per_ndf"])

    out = RESULTS / f"systematics_summary_envelope_{base}.md"
    with out.open("w") as fh:
        fh.write("# Spectral-fit systematic envelope -- 17C(d,p) unbound states\n\n")
        fh.write(f"Binning: **{base}**.  Variants enumerated: "
                 f"{env[0]['n_variants']}.\n\n")

        fh.write("## Inclusive-fit chi^2/NDF grid\n\n")
        fh.write("| cuts \\ contam | nominal | + contam Gauss |\n")
        fh.write("|---|---|---|\n")
        for ct, ct_lab in (("", "loose"), ("strict", "strict")):
            row = f"| **{ct_lab}** |"
            for cm in ("", "contam"):
                v = cells.get((ct, cm))
                row += f" {v:.2f} |" if v is not None else " - |"
            fh.write(row + "\n")
        fh.write("\nThe contaminant Gauss at ~0.45 MeV survives both cut sets, "
                 "dropping chi^2/NDF from ~1.8 (rejected) to ~0.9 (acceptable).\n\n")

        fh.write("## Per-state yield envelope\n\n")
        fh.write("Yield = integrated counts over theta_CM in [10, 40] deg.  "
                 "Envelope spans min..max across all variants.\n\n")
        fh.write("| state | nominal | -% | +% | comment |\n")
        fh.write("|---|---|---|---|---|\n")
        for row in env:
            # frac_lo_pct is signed (negative when env_lo < nominal); same
            # for frac_hi_pct (positive when env_hi > nominal).
            lo = float(row["frac_lo_pct"])      # e.g. -29.9
            hi = float(row["frac_hi_pct"])      # e.g. +0.2
            big = max(abs(lo), abs(hi)) >= 10.0
            mark = "**" if big else ""
            label = STATE_LABELS.get(row["state"], row["state"])
            cmt = ""
            if row["state"] == "ex2.763_1half-":
                cmt = "contam + BW Γ scale"
            elif row["state"] == "ex2.980_3half+":
                cmt = "fit instability between contam and 2.980 BW"
            elif row["state"] == "gs_3half+":
                cmt = "bound state, contam absorbs its high-side tail"
            elif row["state"] == "ex3.661_NEW":
                cmt = "modest, sensitive to contam"
            elif row["state"] == "ex4.231_3half+":
                cmt = "BW Γ scale dominates (narrow BWs grow 4.231)"
            elif row["state"] == "ex4.841_1half+":
                cmt = "BW Γ scale dominates (1.2 MeV state)"
            elif row["state"] == "ex5.91_3half+":
                cmt = "modest"
            elif row["state"] == "ex6.30_5half+":
                cmt = "robust (< 10 %)"
            fh.write(f"| {label} | {mark}{float(row['nominal']):.0f}{mark} "
                     f"| {mark}{lo:+.1f}%{mark} | {mark}{hi:+.1f}%{mark} | {cmt} |\n")
        fh.write("\n**Bold** = envelope width > 10%.\n\n")

        fh.write("## Takeaway for the paper / slides\n\n")
        fh.write("Reading the envelope min..max as the systematic on integrated "
                 "yield, the picture splits the 7 unbound states into three groups.\n\n")
        fh.write("**Large systematic (degeneracy with the contaminant):**\n")
        fh.write("- **2.980 (3/2+)**: -28% / +179%.  The contam Gauss at 0.45 MeV "
                 "overlaps with the low-energy tail of the 2.980 BW; the fit "
                 "cannot fully separate them with the current parametrization. "
                 "Worth flagging as a model-dependence issue, not just a "
                 "statistical uncertainty.\n")
        fh.write("- **2.763 (1/2-)**: -56% / +32%.  Same family of degeneracies "
                 "driven by contam + BW Γ scale interplay.\n\n")
        fh.write("**Moderate systematic (BW Γ scale dominates):**\n")
        fh.write("- **4.231 (3/2+)**: -14% / +48% -- BW Γ = 0.70 (narrow) inflates it.\n")
        fh.write("- **4.841 (1/2+)**: -27% / +11% -- the 1.2 MeV state is the most "
                 "sensitive to the BW Γ scaling because it has the widest range.\n")
        fh.write("- **5.91 (3/2+)**: -7% / +22% -- BW Γ scale.\n\n")
        fh.write("**Robust (statistics dominate):**\n")
        fh.write("- **3.661 NEW**: -15% / +17% -- modest spread despite being a new "
                 "state assignment.\n")
        fh.write("- **6.30 (5/2+)**: -9% / +3% -- robust across all 82 variants.\n")
        fh.write("- **gs 3/2+ (bound)**: -23% / 0% -- contam steals from its "
                 "high-side tail.\n\n")
        fh.write("**Cut variation does NOT remove the contaminant.**  Strict cuts "
                 "(E_ej < 8 MeV + vertex_z in [2, 98] cm) reduce statistics by "
                 "~3x but preserve the 0.45 MeV excess at essentially the same "
                 "centroid; the loose cut set is therefore preferred for the "
                 "final analysis and the contaminant variant is the systematic.\n")
        fh.write("\n**Excluded variants:** the no_3661 family (a state-existence "
                 "test, not a yield systematic) and the orphan contam+relaxsigma "
                 "combination (a known degeneracy between two knobs covering the "
                 "same excess) are removed from the envelope.\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
