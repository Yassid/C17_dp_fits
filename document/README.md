# 16C(d,p)17C analysis report — LaTeX source

Self-contained LaTeX project. Compiles with **pdfLaTeX + BibTeX** (no special
packages beyond a standard TeX Live / Overleaf install).

## Files
- `main.tex` — the full report (~12 pages, 15 figures, 5 tables).
- `refs.bib` — bibliography.
- `figs/` — all figures (our analysis plots + the previous-paper reference
  plots `prev_*`).

## Build (Overleaf)
Create a new Overleaf project and upload this whole `document/` folder, or
drag-and-drop the files. Set the compiler to **pdfLaTeX** and the main document
to `main.tex`. Overleaf runs BibTeX automatically.

## Build (local)
```bash
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Notes
- Editorial placeholders are marked `\note{...}` (renders in red) — search for
  `\note` to find the two spots needing input (physical C²S table, and any
  per-figure refinements).
- Figure provenance: everything except `prev_*.png` was generated in this repo
  (`Codes/`); `prev_inclusive.png` and `prev_10-15.png` are copied verbatim
  from the parent `C16dp_fits` analysis (Movilla & Ayyad) — reference only.
