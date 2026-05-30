#!/usr/bin/env python3
"""Self-contained analysis-overview slide deck (matplotlib PdfPages, 16:9).

No LaTeX toolchain required.  Each slide is a matplotlib page; figure slides
embed the PNGs from plots/ (preserving aspect).  Output:
    Presentations/C17_analysis_overview_auto.pdf
"""
import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

REPO = os.path.join(os.path.dirname(__file__), "..")
PLOTS = os.path.join(REPO, "plots")
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(REPO, "Presentations", "C17_analysis_overview_auto.pdf")

FIG = (13.333, 7.5)
MAROON = "#6f1d46"
ACCENT = "#b3477d"
INK = "#1a1a1a"
GREY = "#555555"

plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK})

LSTATE_LABEL = {
    "ex2.763_1half-": "2.763 1/2-", "ex2.980_3half+": "2.980 3/2+",
    "ex3.661_NEW": "3.661 new", "ex4.231_3half+": "4.231 3/2+",
    "ex4.841_1half+": "4.841 1/2+", "ex5.91_3half+": "5.91 3/2+",
    "ex6.30_5half+": "6.30 5/2+",
}


def read_lscan():
    rows = list(csv.DictReader(open(os.path.join(RESULTS, "fine", "L_scan.csv"))))
    return rows


def slide(title, kicker=None):
    fig = plt.figure(figsize=FIG)
    fig.patch.set_facecolor("white")
    fig.add_artist(Rectangle((0, 0.88), 1, 0.12, color=MAROON, zorder=0))
    fig.add_artist(Rectangle((0, 0.875), 1, 0.006, color=ACCENT, zorder=0))
    fig.text(0.035, 0.927, title, fontsize=23, color="white",
             va="center", fontweight="bold")
    if kicker:
        fig.text(0.965, 0.927, kicker, fontsize=12, color="white",
                 va="center", ha="right", style="italic")
    fig.text(0.035, 0.03,
             r"$^{16}$C(d,p)$^{17}$C  ·  FRIB a1975  ·  Y. Ayyad (USC)",
             fontsize=9, color=GREY)
    return fig


def page(pdf, fig):
    fig.text(0.965, 0.03, str(page.n), fontsize=9, color=GREY, ha="right")
    pdf.savefig(fig)
    plt.close(fig)
    page.n += 1
page.n = 1


def bullets(fig, x, y, items, fs=15, dy=0.062, color=INK):
    for it in items:
        lvl = 0
        s = it
        while s.startswith(">"):
            lvl += 1
            s = s[1:].lstrip()
        marker = "•  " if lvl == 0 else "–  "
        fig.text(x + 0.03 * lvl, y, marker + s, fontsize=fs - 2 * lvl,
                 va="top", ha="left", color=color, wrap=True)
        y -= dy
    return y


def add_image(fig, path, box):
    x, y, w, h = box
    img = plt.imread(path)
    ih, iw = img.shape[0], img.shape[1]
    ar = iw / ih
    fw, fh = fig.get_size_inches()
    bw, bh = w * fw, h * fh
    if bw / bh > ar:
        nw, nh = (bh * ar) / fw, h
        nx, ny = x + (w - nw) / 2, y
    else:
        nw, nh = w, (bw / ar) / fh
        nx, ny = x, y + (h - nh) / 2
    ax = fig.add_axes([nx, ny, nw, nh])
    ax.imshow(img)
    ax.axis("off")


def mono(fig, x, y, rows, fs=12.5, dy=0.045, color=INK):
    for r in rows:
        fig.text(x, y, r, fontsize=fs, va="top", ha="left",
                 family="monospace", color=color)
        y -= dy
    return y


def build():
    L = read_lscan()
    with PdfPages(OUT) as pdf:

        # ---- 1. title -------------------------------------------------
        fig = plt.figure(figsize=FIG); fig.patch.set_facecolor("white")
        fig.add_artist(Rectangle((0, 0), 1, 1, color=MAROON, zorder=0))
        fig.add_artist(Rectangle((0, 0.46), 1, 0.012, color="white", zorder=1))
        fig.text(0.5, 0.66, "Spectroscopy of unbound $^{17}$C states",
                 fontsize=34, color="white", ha="center", fontweight="bold")
        fig.text(0.5, 0.575, "via $^{16}$C(d,p)$^{17}$C at FRIB (a1975)",
                 fontsize=24, color="white", ha="center")
        fig.text(0.5, 0.38, "Angular distributions · absolute normalization · "
                 "$L$ assignments", fontsize=16, color="white", ha="center",
                 style="italic")
        fig.text(0.5, 0.26, "Y. Ayyad   —   IGFAE, Universidade de Santiago de "
                 "Compostela", fontsize=15, color="white", ha="center")
        fig.text(0.5, 0.20, "30 May 2026", fontsize=13, color="white",
                 ha="center")
        fig.text(0.5, 0.08, "spectroscopic reference: Movilla & Ayyad (2026-02)",
                 fontsize=10, color="#e8c9da", ha="center")
        pdf.savefig(fig); plt.close(fig)

        # ---- 2. the experiment ---------------------------------------
        fig = slide("The experiment", "reaction & setup")
        bullets(fig, 0.05, 0.80, [
            r"$^{16}$C beam in the AT-TPC, D$_2$ active target (FRIB a1975).",
            r"188.33 MeV total = 11.77 MeV/u;  $Q_{g.s.}=-1.498$ MeV.",
            r"Transfer populates the bound region + seven unbound resonances.",
            r"Goal: per-state d$\sigma$/d$\Omega$ on an absolute scale",
            r">$\rightarrow$ transferred-$L$ from FRESCO ADWA shape comparison.",
            r"Data reduced with ATTPCROOT (InterpSolver); cont. of C16dp_fits.",
        ], fs=16, dy=0.085)
        fig.text(0.68, 0.80, "Unbound resonances fitted", fontsize=14,
                 fontweight="bold", color=MAROON)
        mono(fig, 0.68, 0.74, [
            "Ex(MeV)   Jpi      Gamma(keV)",
            "-----------------------------",
            "2.763     1/2-          50",
            "2.980     3/2+         200",
            "3.661     new          398",
            "4.231     3/2+         500",
            "4.841     1/2+         300",
            "5.91      3/2+        1500",
            "6.30      5/2+         396",
        ], fs=12.5, dy=0.05)
        page(pdf, fig)

        # ---- 3. pipeline ---------------------------------------------
        fig = slide("Analysis pipeline", "method")
        bullets(fig, 0.05, 0.80, [
            r"1. Spectral fit (C16_pd_AngBins.C): 12 $\times$ 2.5$^\circ$ bins, "
            r"10–40$^\circ$.",
            r">3 Gaussians (bound) + 7 penetrability-modified Breit–Wigners",
            r">+ 1n/2n phase space + linear bg;  analytic Voigt lineshape.",
            r">Inclusive fit pins positions/widths; per-bin frees amplitudes.",
            r"2. Efficiency: 8 Ex-tables, bilinear (Ex,$\theta$); floor+gradient cuts.",
            r"3. Absolute normalization fixed by the $^{16}$C(d,d) elastic channel.",
            r"4. L scan: each AD vs FRESCO L=0–3 ADWA shapes; best L by $\chi^2/\nu$.",
            r"5. Systematics: 82-variant campaign $\rightarrow$ per-state envelope.",
        ], fs=16, dy=0.088)
        page(pdf, fig)

        # ---- 4. inclusive fit ----------------------------------------
        fig = slide("Inclusive $E_x$ spectrum: the global fit", "spectral fit")
        add_image(fig, os.path.join(PLOTS, "plots_inclusive.png"),
                  (0.04, 0.10, 0.58, 0.74))
        bullets(fig, 0.65, 0.80, [
            r"A localized +2–4$\sigma$ excess near 0.4–1.1 MeV is not absorbed "
            r"by the nominal model.",
            r"Adding a contaminant Gaussian (centroid pegs ~0.45 MeV) restores "
            r"an acceptable fit:",
        ], fs=13.5, dy=0.09)
        mono(fig, 0.65, 0.55, [
            "cuts      nominal   +contam",
            "---------------------------",
            "loose      1.78      0.84",
            "strict     1.88      0.97",
            "          (chi2/nu)",
        ], fs=12.5, dy=0.05)
        bullets(fig, 0.65, 0.28, [
            r"Survives strict cuts (0.450$\to$0.451 MeV) $\Rightarrow$ not an "
            r"AT-TPC edge artifact.",
            r"Carried as a model-dependence systematic.",
        ], fs=12.5, dy=0.075)
        page(pdf, fig)

        # ---- 5. bound states -----------------------------------------
        fig = slide("Bound states: angular distributions", "comparison")
        add_image(fig, os.path.join(PLOTS, "bound_states_ad.png"),
                  (0.04, 0.08, 0.60, 0.76))
        bullets(fig, 0.67, 0.80, [
            r"Three bound $^{17}$C states, as a cross-check of the machinery:",
            r">g.s. 3/2$^+$",
            r">0.217 MeV 1/2$^+$",
            r">0.335 MeV 5/2$^+$",
            r"Forward-peaked g.s. and 0.217 shapes consistent with the expected "
            r"$\ell$ content.",
            r"The 5/2$^+$ is only weakly populated — its points sit near the "
            r"detection floor (large bars), shown for completeness.",
            r"Validate efficiency + absolute scale before the unbound region.",
        ], fs=13, dy=0.082)
        page(pdf, fig)

        # ---- 6. normalization idea -----------------------------------
        fig = slide("Absolute normalization: the idea", "normalization")
        fig.text(0.05, 0.80, "Counts $\\rightarrow$ cross section "
                 "(fixed luminosity constant):", fontsize=16, va="top")
        fig.text(0.5, 0.66,
                 r"$\dfrac{d\sigma}{d\Omega}=\dfrac{Y}{N_{beam}\,N_{tgt}\,"
                 r"\Delta\Omega}\cdot 10\;/\;\varepsilon$",
                 fontsize=24, ha="center", va="center")
        bullets(fig, 0.05, 0.52, [
            r"$N_{beam}=1.6146\times10^{5}$,  $N_{tgt}=0.019632$,  "
            r"$\Delta\Omega=2\pi(\cos\theta_{lo}-\cos\theta_{hi})$.",
            r"Cross-check: the $^{16}$C(d,d) elastic channel shares the entrance "
            r"partition (same beam, D$_2$ target, luminosity).",
            r">An absolute OM elastic cross section that reproduces the data "
            r"validates the constant — no free normalization.",
            r"Self-contained partial-wave OM solver (om_elastic.py): Numerov + "
            r"Coulomb matching; vs Rutherford to $10^{-3}$.",
            r">$E_{lab}(d)=23.55$ MeV, $E_{c.m.}=20.92$ MeV.",
        ], fs=15, dy=0.082)
        page(pdf, fig)

        # ---- 7. normalization result ---------------------------------
        fig = slide("Absolute normalization: result", "normalization")
        add_image(fig, os.path.join(PLOTS, "om_elastic_dd.png"),
                  (0.03, 0.08, 0.58, 0.76))
        fig.text(0.64, 0.82, r"Forward-lobe $\langle$data/OM$\rangle$:",
                 fontsize=14, va="top", fontweight="bold")
        mono(fig, 0.64, 0.74, [
            "OMP             chi2/nu   d/OM",
            "------------------------------",
            "ADWA (transfer)   398    0.864",
            "DA1p (global)     105    0.886",
            "free 5-par fit     32    0.978",
        ], fs=12.5, dy=0.05)
        bullets(fig, 0.64, 0.46, [
            r"Absolute scale consistent with unity — no rescale.",
            r"BUT the transfer ADWA potential itself sits ~14% low.",
            r"Honest figure = cross-potential spread $\approx$0.86–0.98 (~12%), "
            r"a GLOBAL scale uncertainty (not in per-point errors).",
        ], fs=12.5, dy=0.082)
        page(pdf, fig)

        # ---- 8. error model ------------------------------------------
        fig = slide("From counts to d$\\sigma$/d$\\Omega$: error model",
                    "uncertainties")
        fig.text(0.5, 0.72,
                 r"$\delta\!\left(\frac{d\sigma}{d\Omega}\right)="
                 r"\frac{d\sigma}{d\Omega}\sqrt{(\delta Y/Y)^2+"
                 r"(\delta\varepsilon/\varepsilon)^2}$",
                 fontsize=23, ha="center", va="center")
        bullets(fig, 0.05, 0.56, [
            r"$\delta Y$ = Minuit amplitude error ($\chi^2$ fit to "
            r"$\sqrt{N}$-weighted histograms) — already carries Poisson stats.",
            r"FIXED (2026-05-30): a previous $\sqrt{1/Y+(\delta Y/Y)^2}$ "
            r"double-counted Poisson ($\sqrt{2}$ inflation).",
            r">e.g. g.s. 10–12.5$^\circ$:  0.433 $\rightarrow$ 0.356 mb/sr.",
            r"Three uncertainty layers, kept separate:",
            r">statistical (above) + systematic (82-variant envelope) "
            r"+ global scale (~12% OMP spread).",
        ], fs=15.5, dy=0.085)
        page(pdf, fig)

        # ---- 9. unbound ADs ------------------------------------------
        fig = slide("Unbound resonances: d$\\sigma$/d$\\Omega$ + "
                    "systematic band", "angular distributions")
        add_image(fig, os.path.join(PLOTS, "unbound_ad_grid.png"),
                  (0.02, 0.10, 0.96, 0.72))
        fig.text(0.5, 0.075, "Black: data ± corrected statistical error.   "
                 "Gray: min–max across the 82 spectral-fit variants.   "
                 r"Fits use $\theta_{c.m.}\leq30^\circ$.",
                 fontsize=11, ha="center", color=GREY)
        page(pdf, fig)

        # ---- 10. L method --------------------------------------------
        fig = slide("$L$ assignment: method", "L values")
        bullets(fig, 0.05, 0.80, [
            r"FRESCO 7-state DWBA with the Johnson–Soper ADWA effective deuteron "
            r"OMP (sum of Koning–Delaroche $p+^{16}$C and $n+^{16}$C at $E_d/2$).",
            r"Four templates L=0,1,2,3 (neutron in $3s_{1/2}$, $2p_{1/2}$, "
            r"$2d_{3/2}/2d_{5/2}$, $1f_{5/2}$); weakly-bound approx $b_e$=0.5 MeV.",
            r"Each AD fit to each L shape (single normalization); best L by "
            r"$\chi^2/\nu$ on $\theta_{c.m.}\leq30^\circ$.",
            r"SHAPE-ONLY: the WBA reproduces the L-dependent shape but not the "
            r"magnitude $\Rightarrow$ $C^2S$ are effective, not publication SFs",
            r">(needs CDCC / pole-residue treatment).",
        ], fs=15.5, dy=0.092)
        page(pdf, fig)

        # ---- 11. FRESCO overlays (NEW) -------------------------------
        fig = slide("$L$ assignment: FRESCO shape overlays", "L values")
        add_image(fig, os.path.join(PLOTS, "plots_fresco_unbound.png"),
                  (0.05, 0.10, 0.90, 0.72))
        fig.text(0.5, 0.075, "Curves = FRESCO ADWA at each state's literature "
                 "single-particle config (normalization-fit, full 10–40$^\\circ$). "
                 "The data-driven all-$L$ scan on $\\theta_{c.m.}\\leq30^\\circ$ "
                 "is the next slide.",
                 fontsize=10.5, ha="center", color=GREY)
        page(pdf, fig)

        # ---- 12. L results (tightened) -------------------------------
        fig = slide("$L$ assignment: results", "L values")
        add_image(fig, os.path.join(PLOTS, "lscan_bars.png"),
                  (0.02, 0.12, 0.55, 0.70))
        fig.text(0.585, 0.82, "best L (corrected errors):", fontsize=12.5,
                 va="top", fontweight="bold")
        mono(fig, 0.585, 0.755, [
            "Ex      best L      next L",
            "-----------------------------",
            "2.763   L=3 (1.41)  L=2 (1.42)",
            "2.980   L=3 (2.73)  L=2 (3.05)",
            "3.661  *L=2 (1.38)  L=3 (1.45)",
            "4.231  *L=2 (0.58)  L=3 (1.38)",
            "4.841  *L=0 (3.09)  L=1 (7.35)",
            "5.91    L=2 (3.55)  L=0 (3.73)",
            "6.30    L=0 (0.34)  L=3 (0.40)",
        ], fs=11, dy=0.044)
        bullets(fig, 0.585, 0.35, [
            r"$\chi^2/\nu$ refreshed after the error fix.",
            r"All best-L picks & tie flags UNCHANGED;",
            r">$\chi^2/\nu$ rose $\times$1.3–1.5 (every L scales together).",
        ], fs=11.5, dy=0.062)
        page(pdf, fig)

        # ---- 13. state-by-state --------------------------------------
        fig = slide("$L$ assignment: state-by-state", "L values")
        bullets(fig, 0.05, 0.80, [
            r"4.231 3/2$^+$ / 4.841 1/2$^+$: decisive (L=2 d-wave; L=0 s-wave) — "
            r"next L worse by factor $\gtrsim$2; match the PDF.",
            r"3.661 (new): L=2 over L=3 (1.38 vs 1.45) — not decisive; possible "
            r"sd–pf intruder.",
            r"2.763 1/2$^-$: three-way near-tie (L=3/2/1); AD too flat to separate.",
            r"2.980 3/2$^+$:  no template fits well; 17.5$^\circ$ AD peak not "
            r"reproduced; fit wants $E_0$ ~120 keV above the SM prior.",
            r"6.30 5/2$^+$:  L=0 has best $\chi^2/\nu$ but is FORBIDDEN "
            r"($0^+\!\otimes s_{1/2}=1/2^+$); allowed L=2/3 are tied.",
            r"5.91 3/2$^+$:  L=2 marginal over L=0; AD poorly constrained "
            r"past 30$^\circ$.",
        ], fs=14.5, dy=0.092)
        page(pdf, fig)

        # ---- 14. systematics -----------------------------------------
        fig = slide("Systematic envelope (82 variants)", "systematics")
        add_image(fig, os.path.join(PLOTS, "yield_envelope_fine.png"),
                  (0.03, 0.10, 0.58, 0.74))
        fig.text(0.64, 0.82, "Integrated-yield envelope:", fontsize=13,
                 va="top", fontweight="bold")
        mono(fig, 0.64, 0.75, [
            "state        -%      +%",
            "----------------------",
            "g.s. 3/2+    -23     +0",
            "2.763        -56    +32",
            "2.980        -28   +179",
            "3.661        -15    +17",
            "4.231        -14    +48",
            "4.841        -27    +11",
            "5.91          -8    +22",
            "6.30          -9     +3",
        ], fs=12, dy=0.045)
        bullets(fig, 0.64, 0.30, [
            r"Largest model-dependence: 2.763 & 2.980 (contaminant degeneracy).",
            r"6.30 robust (<10%).",
        ], fs=12, dy=0.07)
        page(pdf, fig)

        # ---- 15. summary ---------------------------------------------
        fig = slide("Summary", "")
        bullets(fig, 0.05, 0.80, [
            r"Per-state absolute d$\sigma$/d$\Omega$ for 3 bound + 7 unbound "
            r"$^{17}$C states from $^{16}$C(d,p) at 11.77 MeV/u.",
            r"Absolute scale validated vs $^{16}$C(d,d) elastic; normalization "
            r"uncertainty = ~12% OMP spread (global).",
            r"Error model corrected: per-point stats no longer double-counted; "
            r"statistical + systematic + scale layers separated.",
            r"L assignments robust to the fix: 4.231 (L=2) & 4.841 (L=0) "
            r"decisive; 3.661/2.763/5.91 ambiguous; 2.980 & 6.30 problematic.",
            r"Deliverable (this repo): yields, d$\sigma$/d$\Omega$, systematic "
            r"envelope. Final L / $C^2S$ done downstream (CDCC-grade).",
        ], fs=15.5, dy=0.10)
        page(pdf, fig)

        # ---- 16. open items ------------------------------------------
        fig = slide("Open items", "")
        bullets(fig, 0.05, 0.80, [
            r"Fold the ~12% elastic-normalization scale into the published "
            r"d$\sigma$/d$\Omega$ / yield bands (common factor).",
            r"Continuum-form tests (quadratic/exp bg, PS-shape) to show the "
            r"0.45 MeV contaminant is not bg flexibility.",
            r"Physically identify the 0.45 MeV excess (vertex_z, "
            r"$\theta_{lab}$ fingerprint of a non-D$_2$ channel?).",
            r"CDCC / pole-residue treatment for publication-grade absolute "
            r"$C^2S$.",
        ], fs=16, dy=0.10)
        fig.text(0.05, 0.16, "Refs: Pereira-Lopez et al., PLB 811 (2020) 135939; "
                 "Movilla & Ayyad (2026); Pang et al., arXiv:1606.01507 (DA1p).",
                 fontsize=10, color=GREY)
        page(pdf, fig)

        # ---- 17. backup: full chi2/nu grid (NEW) ---------------------
        fig = slide("Backup: full $\\chi^2/\\nu$ grid", "L values")
        hdr = "state         L=0     L=1     L=2     L=3    best"
        lines = [hdr, "-" * len(hdr)]
        for r in L:
            lab = LSTATE_LABEL.get(r["state"], r["state"])
            vals = [float(r[f"c2v_L{i}"]) for i in range(4)]
            best = int(r["bestL"])
            cells = []
            for i, v in enumerate(vals):
                cells.append(("*" if i == best else " ") + f"{v:5.2f}")
            tie = ("  (tie L=%s)" % r["tieL"]) if r["tie"] == "1" else ""
            lines.append(f"{lab:<11s}" + " ".join(cells) +
                         f"   L={best}{tie}")
        mono(fig, 0.06, 0.78, lines, fs=14, dy=0.058)
        fig.text(0.06, 0.20,
                 "* = best L.  Fits on theta_cm <= 30 deg, corrected stat errors "
                 "(2026-05-30).\nA tie = runner-up within 0.3 in chi2/nu "
                 "(shapes indistinguishable at current statistics).",
                 fontsize=11, color=GREY, va="top")
        page(pdf, fig)

        # ---- 18. backup: non-resonant background --------------------
        fig = slide("Backup: non-resonant background treatment", "background")
        add_image(fig, os.path.join(PLOTS, "fine", "plots_bins.png"),
                  (0.015, 0.10, 0.60, 0.74))
        add_image(fig, os.path.join(PLOTS, "bg_vs_angle.png"),
                  (0.63, 0.46, 0.36, 0.38))
        bullets(fig, 0.63, 0.40, [
            r"In the fit: 1n PS + linear bg + 2n PS, amplitudes free in EVERY "
            r"angular bin (PS templates rebuilt per bin).",
            r"Left: per-bin $E_x$ spectra — peak-dominated forward; from "
            r"~30$^\circ$ the transfer peaks fade and the continuum dominates.",
            r"Right: scaling each term 0.3$\times$–3$\times$ moves "
            r"d$\sigma$/d$\Omega$ by <2% (bg), ~0% (PS) — flat in angle, no "
            r"rise past 30$^\circ$.",
            r"Also in the 82-variant envelope; L-scan uses "
            r"$\theta\leq30^\circ$. Back-angle limit is stats/efficiency, not bg.",
        ], fs=11, dy=0.058)
        page(pdf, fig)

    print("wrote", OUT, "(", page.n - 1, "content slides + title )")


if __name__ == "__main__":
    build()
