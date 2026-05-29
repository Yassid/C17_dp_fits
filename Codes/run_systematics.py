"""
Driver for spectral-fit systematics on the 17C unbound-states analysis.

For each variant in VARIANTS, sets env vars, runs C16_pd_AngBins.C, runs
best_L_scan.py, then reads back results/<variant>/L_scan.csv and writes
a single summary at results/systematics_summary_<binning>.{csv,md}.

Run from Codes/:
    python3 run_systematics.py             # fine binning, all variants
    python3 run_systematics.py --only nominal,no_3661,bw_g1.3
    python3 run_systematics.py --binning coarse2 --only nominal,bw_g0.7
    python3 run_systematics.py --dry        # just print what would run
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
RESULTS = REPO / "results"

# Scheme: BINNING name (passed to best_L_scan.py); macro arg 0/1/2.
BIN_TO_ARG = {"fine": 0, "coarse2": 1, "coarse4": 2}

# Variant = (name, env-var dict).  env vars not listed default to nominal.
# Each variant lands in results/<binning>_<name>/ (nominal goes to <binning>/).
VARIANTS = [
    # nominal must be first.
    ("nominal", {}),

    # ---- BW-width scaling: all 7 unbound widths +/- 30% and +/- 15% ----
    ("bw_g0.70", {"BW_GAMMA_SCALE": "0.70"}),
    ("bw_g0.85", {"BW_GAMMA_SCALE": "0.85"}),
    ("bw_g1.15", {"BW_GAMMA_SCALE": "1.15"}),
    ("bw_g1.30", {"BW_GAMMA_SCALE": "1.30"}),

    # ---- State counting: drop the 3.661 NEW state and refit ----
    ("no_3661", {"STATE_DROP": "0,0,1,0,0,0,0"}),

    # ---- Phase-space / linear-bg scales: each amplitude floor x0.3 / x3.0 ----
    ("ps1n_x0.30", {"PS1N_SCALE": "0.30"}),
    ("ps1n_x3.00", {"PS1N_SCALE": "3.00"}),
    ("ps2n_x0.30", {"PS2N_SCALE": "0.30"}),
    ("ps2n_x3.00", {"PS2N_SCALE": "3.00"}),
    ("bg_x0.30",   {"BG_LINEAR_SCALE": "0.30"}),
    ("bg_x3.00",   {"BG_LINEAR_SCALE": "3.00"}),

    # ---- Per-state penetrability-L overrides ----
    # Each tweaks BW_L for one state only; everything else stays nominal.
    # nominal kBW_L = {1,2,2,2,2,2,2}.  Alternatives tried:
    #   state 0 (2.763 1/2-): try L=3
    #   state 1 (2.980 3/2+): try L=0 and L=3
    #   state 2 (3.661 NEW):  try L=3 (data preference) and L=0
    #   state 3 (4.231 3/2+): try L=0
    #   state 4 (4.841 1/2+): try L=0 (decisive L=0 in data)
    #   state 5 (5.91 3/2+):  try L=3
    #   state 6 (6.30 5/2+):  try L=3 (PDF candidate)
    ("L_s0_L3", {"BW_L_OVERRIDE": "3,2,2,2,2,2,2"}),
    ("L_s1_L0", {"BW_L_OVERRIDE": "1,0,2,2,2,2,2"}),
    ("L_s1_L3", {"BW_L_OVERRIDE": "1,3,2,2,2,2,2"}),
    ("L_s2_L0", {"BW_L_OVERRIDE": "1,2,0,2,2,2,2"}),
    ("L_s2_L3", {"BW_L_OVERRIDE": "1,2,3,2,2,2,2"}),
    ("L_s3_L0", {"BW_L_OVERRIDE": "1,2,2,0,2,2,2"}),
    ("L_s4_L0", {"BW_L_OVERRIDE": "1,2,2,2,0,2,2"}),
    ("L_s5_L3", {"BW_L_OVERRIDE": "1,2,2,2,2,3,2"}),
    ("L_s6_L3", {"BW_L_OVERRIDE": "1,2,2,2,2,2,3"}),
]

# Axis variants -- multiplied against the base list when --axes is passed.
# Each entry is (tag, env-dict).  Empty tag = no name change (default cell).
CUT_AXIS = [
    ("",       {}),                              # nominal: polar>=90 only
    ("strict", {"CUT_SET": "strict"}),           # adds E_ej<8 + vertex_z in [2,98] cm
]
CONTAM_AXIS = [
    ("",       {}),                              # no contaminant Gauss
    ("contam", {"CONTAMINANT_ENABLE": "1"}),     # 4th Gauss in [0.45, 1.20] MeV
]


def build_cross_product(base_variants, axes: list[str]) -> list[tuple[str, dict]]:
    """Expand each base variant across the requested axes.

    `axes` is a subset of {"cuts", "contam"}.  When empty the returned list
    is exactly base_variants.  Default-cell (empty tag) cells of all chosen
    axes collapse onto the existing variant names so the 21 nominal-cell
    results are not re-run with a different name.
    """
    if not axes:
        return list(base_variants)
    chosen_axes = []
    if "cuts" in axes:   chosen_axes.append(CUT_AXIS)
    if "contam" in axes: chosen_axes.append(CONTAM_AXIS)
    out = []

    def rec(combo, env_acc, tags):
        if not combo:
            for base_name, base_env in base_variants:
                env = {**env_acc, **base_env}
                parts = [t for t in tags if t]
                if base_name != "nominal":
                    parts.append(base_name)
                name = "_".join(parts) if parts else "nominal"
                out.append((name, env))
            return
        for tag, axis_env in combo[0]:
            rec(combo[1:], {**env_acc, **axis_env}, tags + [tag])

    rec(chosen_axes, {}, [])
    seen, dedup = set(), []
    for n, e in out:
        if n in seen: continue
        seen.add(n); dedup.append((n, e))
    return dedup


STATE_KEYS = [
    "ex2.763_1half-",
    "ex2.980_3half+",
    "ex3.661_NEW",
    "ex4.231_3half+",
    "ex4.841_1half+",
    "ex5.91_3half+",
    "ex6.30_5half+",
]


def variant_dir(binning: str, vname: str) -> Path:
    sub = binning if vname == "nominal" else f"{binning}_{vname}"
    return RESULTS / sub


def run_macro(binning: str, vname: str, env_extra: dict, log: Path) -> bool:
    macro_arg = BIN_TO_ARG[binning]
    env = os.environ.copy()
    env["VARIANT_NAME"] = vname
    for k, v in env_extra.items():
        env[k] = v
    cmd = ["root", "-l", "-b", "-q", f"C16_pd_AngBins.C++({macro_arg})"]
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as fh:
        fh.write(f"# variant={vname}  env_extra={env_extra}\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT)
    return rc == 0


def run_lscan(binning: str, vname: str, log: Path) -> bool:
    env = os.environ.copy()
    env["BINNING"] = binning if vname == "nominal" else f"{binning}_{vname}"
    cmd = ["python3", "best_L_scan.py"]
    with log.open("a") as fh:
        fh.write("\n# --- best_L_scan ---\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT)
    return rc == 0


def read_lscan(binning: str, vname: str) -> dict:
    """Return {state: dict(...)} from results/<dir>/L_scan.csv."""
    csv_path = variant_dir(binning, vname) / "L_scan.csv"
    out = {}
    if not csv_path.exists():
        return out
    with csv_path.open() as fh:
        for row in csv.DictReader(fh):
            out[row["state"]] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", default="fine", choices=list(BIN_TO_ARG))
    ap.add_argument("--only", default="", help="comma list of variant names; default = all")
    ap.add_argument("--axes", default="",
                    help="comma list subset of {cuts,contam}; multiplies the variant set")
    ap.add_argument("--dry", action="store_true", help="just print plan")
    args = ap.parse_args()

    axes = [a for a in args.axes.split(",") if a]
    all_variants = build_cross_product(VARIANTS, axes)
    selected = [v for v in all_variants
                if not args.only or v[0] in args.only.split(",")]
    if not selected:
        print(f"no variants matched --only {args.only!r}", file=sys.stderr)
        sys.exit(2)

    print(f"binning = {args.binning}   variants = {len(selected)}")
    for name, env in selected:
        print(f"  {name:<14}  {env}")
    if args.dry:
        return

    logs_dir = HERE / "syst_logs" / args.binning
    logs_dir.mkdir(parents=True, exist_ok=True)
    timings = []
    for name, env in selected:
        log = logs_dir / f"{name}.log"
        t0 = time.time()
        print(f"\n[ {name:<14} ]  macro ...", flush=True)
        ok = run_macro(args.binning, name, env, log)
        if not ok:
            print(f"  macro failed -- see {log}")
            continue
        print(f"  L-scan ...", flush=True)
        ok = run_lscan(args.binning, name, log)
        if not ok:
            print(f"  L-scan failed -- see {log}")
            continue
        dt = time.time() - t0
        timings.append((name, dt))
        print(f"  done in {dt:.1f} s")

    # Aggregate.
    print("\n--- aggregating ---")
    rows = []
    nom = read_lscan(args.binning, "nominal")
    if not nom:
        print("nominal L_scan.csv missing -- cannot diff against nominal")
    for name, _env in selected:
        data = read_lscan(args.binning, name)
        for s in STATE_KEYS:
            r = data.get(s, {})
            n = nom.get(s, {})
            best_L  = r.get("bestL", "")
            best_n  = n.get("bestL", "")
            tie     = r.get("tie", "0")
            tieL    = r.get("tieL", "-1")
            row = {
                "variant": name,
                "state":   s,
                "bestL":   best_L,
                "tie":     tie,
                "tieL":    tieL,
                "flipped_vs_nominal": ("1" if (best_L != best_n and best_n) else "0"),
                "c2v_L0":  r.get("c2v_L0", "nan"),
                "c2v_L1":  r.get("c2v_L1", "nan"),
                "c2v_L2":  r.get("c2v_L2", "nan"),
                "c2v_L3":  r.get("c2v_L3", "nan"),
            }
            rows.append(row)

    out_csv = RESULTS / f"systematics_summary_{args.binning}.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}")

    # Markdown table -- per state, list each variant's best L and chi2/nu.
    out_md = RESULTS / f"systematics_summary_{args.binning}.md"
    L_LABEL = ["L=0", "L=1", "L=2", "L=3"]
    with out_md.open("w") as fh:
        fh.write(f"# Spectral-fit systematics ({args.binning})\n\n")
        fh.write("Best L per state under each variant (* = flipped vs nominal). "
                 "chi2/nu shown for the variant's best L.\n\n")
        # Header
        names = [n for n, _ in selected]
        fh.write("| state | " + " | ".join(names) + " |\n")
        fh.write("|" + "---|" * (len(names) + 1) + "\n")
        # Per state row
        by_state = {s: {r["variant"]: r for r in rows if r["state"] == s}
                    for s in STATE_KEYS}
        for s in STATE_KEYS:
            cells = []
            for n in names:
                r = by_state[s].get(n, {})
                bL = r.get("bestL", "")
                if bL in ("", "-1"):
                    cells.append("none")
                    continue
                bLi = int(bL)
                c2v = r.get(f"c2v_L{bLi}", "nan")
                try:
                    c2v_f = float(c2v)
                    cell = f"{L_LABEL[bLi]} ({c2v_f:.2f})"
                except ValueError:
                    cell = L_LABEL[bLi]
                if r.get("tie") == "1":
                    cell += "*tie"
                if r.get("flipped_vs_nominal") == "1":
                    cell = "**" + cell + "**"
                cells.append(cell)
            fh.write("| " + s + " | " + " | ".join(cells) + " |\n")
        # Timing footer
        fh.write("\n## Run times\n\n")
        for n, dt in timings:
            fh.write(f"- {n}: {dt:.1f} s\n")
    print(f"wrote {out_md}")


if __name__ == "__main__":
    main()
