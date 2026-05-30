# TODO — ¹⁶C(d,p)¹⁷C analysis (open items)

Snapshot 2026-05-30, after the hardening + deck work (branch
`harden-and-deck-2026-05-30`). Done this round: per-point stat-error fix +
CSV reconciliation, absolute-normalization re-sync + ±12% correlated scale
band (`results/{dsdo,yield}_bands_fine.csv`), L-scan refresh, 18-slide
overview deck, and the non-resonant-background diagnostic. The threads below
are what's left.

## Physics

1. **2.980 MeV 3/2⁺ — no template fits.**
   Best χ²/ν is only 2.73 (L=3); the 17.5° peak in the AD is not reproduced,
   and the spectral fit wants E₀ ≈ +120 keV above the shell-model prior
   (see [[feedback-shellmodel-priors]]). Check: two unresolved states?
   centroid/width misassignment? contaminant leakage? Decide the final quoted
   L (or flag "ambiguous").

2. **6.30 MeV 5/2⁺ — selection-rule conflict.**
   Lowest χ²/ν is L=0 (0.34) but 0⁺⊗s₁/₂ = 1/2⁺ forbids L=0 for a 5/2⁺.
   Quote the allowed **L=2 (0.41) / L=3 (0.40)** tie; do not report L=0.

3. **Ambiguous L: 3.661 (new), 2.763, 5.91.**
   Document as ambiguous. Final L is assigned downstream by the theorist
   (see [[feedback-l-assignment-defer-to-theory]]); this repo delivers
   yields + dσ/dΩ + systematic envelope, not the final L.

4. **0.45 MeV contaminant — physical ID (it's a BACK-ANGLE feature).**
   Already handled *in the analysis* (4th Gaussian, free amplitude per bin,
   carried as the `contam` systematic; survives strict cuts → not an AT-TPC
   edge artifact). **Its angular profile is already telling**: from
   `results/fine_contam/yields.csv`, the contaminant is pinned at the fit floor
   (~0) for θ_cm ≲ 30° and carries essentially all its strength at the back —
   ≈34, 42, 91, 9 counts at 31/34/36/39° vs ~0 forward (**~87% of the total at
   θ_cm ≥ 31°**; one poorly-constrained 22-count blip at 16°). A feature at a
   *fixed apparent* Eₓ ≈ 0.45 MeV that only appears at back angles points to a
   kinematic mis-ID / different reaction locus, not a real ¹⁷C level.
   **Next:** vertex_z and θ_lab fingerprint of the 0.4–1.1 MeV window to confirm
   a non-D₂ channel (¹²C(d,p) backing, beam contaminant); if clean, convert the
   systematic into a cut.

5. **Publication-grade absolute C²S.**
   Current C²S are *shape-only* (weakly-bound approximation, bₑ=0.5 MeV):
   correct L-dependent angular shape, artificial magnitude. A continuum-binned
   (CDCC-like) or pole-residue treatment is needed for real spectroscopic
   factors.

## Deliverables / hygiene

6. **Normalization scale (~12%).**
   Folded as a *correlated* outer band in `dsdo_bands_fine.csv` /
   `yield_bands_fine.csv` (`fold_scale_band.py`). Confirm this is how it should
   appear in the paper figures (outer band, not per-point) and propagate to the
   published AD plots.

7. **Continuum-form tests.**
   Quadratic / exponential bg and phase-space-shape variants, to show the
   0.45 MeV excess is not just background flexibility. Variants planned, not run.

## Provenance notes

- L-scan χ²/ν are post-error-fix (rose ×1.3–1.5 vs the old inflated bars);
  best-L picks and tie flags unchanged.
- Elastic data `~/Downloads/16C_dd_gs.txt` dated 2026-05-29; `om_elastic.py`
  numbers reproduce live (ADWA 0.864 / DA1p 0.886 / fit 0.978).
- Deck is regenerable: `deck_figs.py` → `make_deck.py` (PDF) / `make_pptx.py`
  (PPTX); Beamer source `Presentations/C17_analysis_overview.tex`.
