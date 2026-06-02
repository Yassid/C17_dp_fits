# 17C unbound-states analysis from 16C(d,p)

End-to-end pipeline for extracting angular distributions of the 7 unbound
resonances of 17C populated in the 16C(d,p)17C reaction at 11.77 MeV/u
(experiment a1975 at FRIB / AT-TPC), and fitting them with FRESCO ADWA
single-particle shapes to assign orbital angular momentum L.

Continuation of the analysis in
[github.com/Yassid/C16dp_fits](https://github.com/Yassid/C16dp_fits), with the
PDF [`Presentations/Analysis_16C_dp.pdf`](Presentations/Analysis_16C_dp.pdf)
(Movilla & Ayyad, 2026-02) as the spectroscopic reference.

## Reaction and beam

- 16C beam: 188.33 MeV total = **11.77 MeV/u** (Q(g.s.) = −1.498 MeV)
- D2 target in the AT-TPC
- 7 unbound 17C states fitted as penetrability-modified Breit-Wigners:
  | Ex (MeV) | J^π          | Γ (keV)   |
  |----------|--------------|-----------|
  | 2.763    | 1/2−         | 50        |
  | 2.980    | 3/2+         | 200       |
  | 3.661    | NEW          | 398       |
  | 4.231    | 3/2+         | 500       |
  | 4.841    | 1/2+         | 300       |
  | 5.91     | 3/2+         | 1500      |
  | 6.30     | 5/2+         | 396       |

## Pipeline

1. `Codes/C16_pd_AngBins.C`
   Reads the InterpSolver ROOT files, computes Q-corrected Ex per event,
   fills 12 × 2.5° angular bins (10–40°) plus an inclusive 10–40° histogram.
   Fits the inclusive Ex spectrum with 3 Gaussians (bound states) + 7
   penetrability-modified Breit-Wigners + 1n & 2n phase space + linear bg,
   then refits each angular bin with positions / widths frozen and only
   amplitudes free. Writes per-state integrated counts (`results/yields.csv`)
   and dσ/dΩ in mb/sr (`results/dsdo.csv`). Plot panels in `plots/`.

   BW lineshape uses `BWModificada` (energy-dependent Γ via tabulated
   penetrabilities in [`Penetrabilities_def/`](Penetrabilities_def/) for
   L=0..3) and an **analytic Voigt** convolution
   (`Amp · 2π · TMath::Voigt(E−E_eff, σ, Γ_eff, 4)`), ~1000× faster than the
   1000-point numerical integral in the original macro.

2. `Codes/plot_angular_distros.py`
   Reads `results/dsdo.csv` and plots dσ/dΩ vs θ_CM for each of the 7 unbound
   resonances (plus the 3 bound Gaussians, kept for cross-check).

3. `fresco/inputs/dp17C_unbound_adwa.nin`
   FRESCO 7-state DWBA with the Johnson-Soper ADWA effective deuteron OMP
   (sum of Koning-Delaroche p+16C and n+16C at E_d/2). All 7 states treated
   in the weakly-bound approximation with be = 0.5 MeV and the n
   configuration appropriate to the literature J^π
   (2p1/2, 2d3/2, 3s1/2, 2d5/2). Output curves in `fresco/outputs/`.

4. `Codes/fresco_overlay_unbound.py`
   Overlays the FRESCO shape for each state against the extracted AD and
   fits a single normalisation (effective C²S).

5. `fresco/inputs/gen_L_templates.py`
   Builds four FRESCO inputs (`dp17C_L{0,1,2,3}.nin`) — one per orbital
   angular momentum of the transferred neutron, with all 7 Ex values at the
   matching (n, l, j). With `--run`, executes them via `run_fresco.sh`.

6. `Codes/best_L_scan.py`
   For each of the 7 states fits the extracted AD against the L=0,1,2,3
   FRESCO shapes and picks the best L by χ²/ν. Produces the
   `plots/plots_L_scan.png` 2 × 4 summary panel.

## Headline results

L-scan, fit on θ_CM ≤ 30° (kinematic-edge bins excluded), weakly-bound
approximation (be = 0.5 MeV in the WS for the n+16C form factor).
Numbers below use the post-May-20 efficiency baseline (inclusive
spectral-fit χ²/ν = 1.78, from `results/fine/L_scan.csv`).

*χ²/ν refreshed 2026-05-30 after the dσ/dΩ stat-error fix (see *Absolute
normalization* below): all seven L picks are unchanged and no near-tie
flipped; the absolute χ²/ν rose ×1.3–1.5 because every L scales by the
same factor.*

| Ex (MeV) | Lit J^π | best L (χ²/ν) | next L (χ²/ν)          | comment                                                            |
|----------|---------|---------------|------------------------|--------------------------------------------------------------------|
| 2.763    | 1/2−    | L=3 (1.41)    | L=2 (1.42), L=1 (1.61) | three-way near-tie; AD too flat to distinguish                     |
| 2.980    | 3/2+    | L=3 (2.73)    | L=2 (3.05)             | all L fits poor; 17.5° peak still not reproduced                   |
| 3.661    | NEW     | **L=2 (1.38)**| L=3 (1.45)             | L=2/L=3 ambiguous; possible sd-pf intruder                         |
| 4.231    | 3/2+    | **L=2 (0.58)**| L=3 (1.38)             | decisive d-wave; matches PDF                                       |
| 4.841    | 1/2+    | **L=0 (3.09)**| L=1 (7.35)             | decisive s1/2 (L=1/2/3 all factor ≥2 worse)                        |
| 5.91     | 3/2+    | L=2 (3.55)    | L=0 (3.73)             | L=2/L=0 near-tie; AD shape poorly constrained                      |
| 6.30     | 5/2+    | L=0 (0.34)    | L=3 (0.40), L=2 (0.41) | L=0 unphysical for 5/2⁺ (s1/2⊗0⁺ = 1/2⁺); L=2/3 essentially tied   |

The absolute C²S_eff values are **shape-only**: the weakly-bound
approximation gives the correct L-dependent angular shape but the magnitude
is artificial (FRESCO sp-strength=1 peaks at ~10²–10³ mb/sr in this
prescription). A continuum-binned (CDCC-like) or pole-residue treatment is
needed for publication-quality absolute SFs.

## Absolute normalization and dσ/dΩ uncertainties

**Counts → cross section** (`Codes/C16_pd_AngBins.C`):
σ(θ) = Y / (N_beam · N_target · ΔΩ) · 10 / ε, with N_beam = 161460.66,
N_target = 0.019632 (the calibrated beam×target product distilled from the
bound-state cell-7 absolute-norm calc), ΔΩ = 2π(cosθ_lo − cosθ_hi), and ε the
(Ex, θ)-interpolated efficiency. Verified: the g.s. 10–12.5° point reproduces
exactly from `yields.csv` (Y = 11.55 → σ = 0.829 mb/sr).

**Absolute-scale check** (`Codes/om_elastic.py`, re-run 2026-05-30 against
`~/Downloads/16C_dd_gs.txt`): the d+16C elastic forward lobe shares the
entrance partition (same beam, target, luminosity) as the transfer, so an
absolute OM elastic cross section that matches the data fixes the
(N_beam · N_target) constant. Error-weighted ⟨data/OM⟩ over the forward lobe:

| OMP                        | χ²/ndf | ⟨data/OM⟩ |
|----------------------------|--------|-----------|
| ADWA (the transfer OMP)    | 398    | 0.864     |
| DA1p (light-target global) | 105    | 0.886     |
| free 5-param fit           | 32     | 0.978     |

The free-fit potential reproduces the absolute scale to ~2% (⟨data/OM⟩ =
0.978 ± 0.043), but the **ADWA potential actually used in the transfer DWBA
gives 0.864** — the data sit ~14 % below it. The honest normalization
uncertainty is the **cross-potential spread, ≈0.86–0.98 (~12 %)**, to be
quoted as a global scale uncertainty; it is NOT folded into the per-point
`dsdo_err`. (Earlier notes citing "ADWA ≈ 0.98, χ²/ndf ≈ 1.8" predate the
current `16C_dd_gs.txt`, dated 2026-05-29, and no longer reproduce.)

**Per-point errors** (`dsdo_err_mbsr`): statistical only,
σ·√((dY/Y)² + (ε_err/ε)²), with dY the Minuit amplitude error (a χ²-fit to
√N-weighted histograms, so it already carries the Poisson statistics). A
previous σ_raw·√(1/Y + (dY/Y)²) **double-counted the statistics** (√2
inflation in the Poisson limit) — fixed in the macro 2026-05-30; the already-
generated CSVs were reconciled in place by `Codes/fix_stat_err.py` (9714 rows
across 90 files) and `Codes/fix_band_stat_err.py` (the systematic band)
without re-running ROOT. Per-state **systematics** are the 82-variant min/max
envelope in `results/systematics_dsdo_band_fine.csv`. The ~12% absolute-scale
uncertainty is **correlated** across all states/angles (it shifts the absolute
scale, not the relative picture), so it is carried as a separate outer band
rather than added per-point; combined bands are tabulated in
`results/dsdo_bands_fine.csv` and `results/yield_bands_fine.csv`
(`Codes/fold_scale_band.py`).

## Running it

ROOT environment from FairSoft (assumes the user's local install):

```bash
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
```

Spectral-fitting macro (data path is `/home/yassid/C16_dp/C16_dp/InterpSolver_root/`,
edit `dataDir` in `C16_pd_AngBins.C` if your data is elsewhere):

```bash
cd Codes
root -l -b -q 'C16_pd_AngBins.C++'        # → results/{yields,dsdo}.csv, plots/plots_{inclusive,bins}.png
python plot_angular_distros.py            # → plots/plots_angular_distros.png
```

FRESCO calculations (assumes `fresco` is in `$PATH`):

```bash
cd fresco
./run_fresco.sh dp17C_unbound_adwa
python inputs/gen_L_templates.py --run    # runs all four dp17C_L{0..3}.nin
```

Shape analysis:

```bash
cd Codes
python fresco_overlay_unbound.py          # → plots/plots_fresco_unbound.png
python best_L_scan.py                     # → plots/plots_L_scan.png
```

## Repository layout

```
Codes/                  ROOT and Python analysis
  C16_pd_AngBins.C      spectral fits → yields.csv, dsdo.csv
  plot_angular_distros.py
  fresco_overlay_unbound.py
  best_L_scan.py
Penetrabilities_def/    tabulated n+16C penetrabilities, L = 0..3
Phase_Space/            1n and 2n phase-space ROOT templates
Presentations/          Analysis_16C_dp.pdf (Movilla/Ayyad reference)
fresco/
  run_fresco.sh         elastic+transfer extractor (fixed: keeps both
                        deuteron-elastic state #1 and proton-transfer
                        state #1 separate)
  inputs/
    kd_adwa.py          builds Johnson-Soper ADWA OMP from Koning-Delaroche
    gen_L_templates.py  emits dp17C_L{0..3}.nin
    dp17C_unbound_adwa.nin
    dp17C_L{0,1,2,3}.nin
  outputs/              per-state per-L dσ/dΩ curves
results/                yields.csv, dsdo.csv
plots/                  PNG output of all analysis steps
```

## Caveats and open items

- **Stats vs. published PDF**: closely matched once the run list and cuts
  align with the parent macro `C16_pd_ana.C` (45 runs from 0016 to 0098,
  plus 0102/0103; cut set: `polar ≥ 90°` only — the `_ang_dist.C`
  version's `E_ej < 8 MeV` and `vertex_z ∈ [2,98] cm` cuts drop ~half the
  events and were not used for the PDF presentation).
- **35–40° bins**: bound-state Gaussian amplitudes hit the fit floor at the
  AT-TPC kinematic edge — drop these bins when comparing bound-state ADs.
- **2.980 3/2+**: no L template fits the data well (best χ²/ν = 2.73 at
  L=3, but expected L=2 only marginally worse at 3.05); the 17.5° peak in
  the AD is not reproduced. Likely background contamination under the BW
  or a centroid/width misassignment — the spectral fit also wants E₀
  ~120 keV above the shell-model value (see [[feedback-shellmodel-priors]]).
- **3.661 NEW**: L=2 narrowly preferred over L=3 (1.38 vs 1.45) with the
  new baseline — flipped from pre-efficiency-fix where L=3 led. Neither
  decisive at current statistics.
- **4.841 1/2+**: L=0 is still the decisive pick (next L worse by factor
  ≥2 in χ²/ν), but the absolute χ²/ν climbed from 0.30 to 3.09 after the
  efficiency correction — worth a residual check in the AD, since the
  spectral fit also resolves a Γ ~3× the shell-model value here.
- **5.91 3/2+ / 6.30 5/2+**: kinematic-edge bins drop most of the signal
  past 30°, so the AD is fit on only ~8 effective points. 6.30 picks
  L=0 by χ²/ν but that is forbidden by 0⁺⊗s1/2 = 1/2⁺; the L=2 d-wave
  fit (χ²/ν = 0.41) is essentially tied and physically allowed.
- **Absolute SFs**: not extracted here. The shape-only weakly-bound
  approximation is sufficient for L assignment, not for spectroscopic-factor
  publication.

## Acknowledgements / references

- ATTPCROOT (a1975 reduction → `parquettree`) — input data files at
  `/home/yassid/C16_dp/C16_dp/InterpSolver_root/`.
- Daniel Movilla & Yassid Ayyad, *16C(d,p) analysis*, FRIB internal slide
  deck, 2026-02-17.
- Pereira-Lopez et al., **Phys. Lett. B 811 (2020) 135939** — paper-energy
  17C(d,p) reference data and benchmark spectroscopic factors.
