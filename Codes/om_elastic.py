#!/usr/bin/env python3
"""
Self-contained partial-wave optical-model elastic-scattering solver.

Purpose: validate the absolute normalization of the 16C(d,p)17C analysis using
the measured d+16C elastic angular distribution (16C_dd_gs.txt).  The elastic
channel shares the entrance partition -- same beam, same D2 target, same
luminosity -- as the transfer, so an absolute optical-model elastic cross
section that reproduces the data confirms the (Nbeam * Ntarget) constant baked
into C16_pd_AngBins.C.

Method (textbook, e.g. Thompson & Nunes "Nuclear Reactions for Astrophysics"):
  * complex optical potential U(r) = V_C - V f_V - i(W_V f_W + W_D g_D)
  * integrate the radial Schrodinger equation per partial wave L (Numerov)
  * match to Coulomb functions F_L,G_L at r_match -> nuclear S-matrix S_L
  * f(theta) = f_C(theta) + (1/2ik) sum_L (2L+1) e^{2 i sigma_L}(S_L-1) P_L
  * dsigma/dOmega = |f|^2  [fm^2/sr] * 10 -> mb/sr

Default configuration: d+16C, Elab(d)=23.55 MeV, Johnson-Soper ADWA potential
(kp=2 in fresco/inputs/dp17C_L0.nin).  No spin-orbit (ADWA central-only), so
the deuteron spin is not resolved -- adequate for an absolute-scale check.

Coulomb functions via mpmath; everything else numpy.
"""
from __future__ import annotations

import argparse
import math
import cmath

import numpy as np
import mpmath as mp
from scipy.special import eval_legendre

# ---- physical constants -------------------------------------------------
E2 = 1.43996_4548          # e^2 in MeV*fm (Coulomb constant)
HBAR2_2U = 20.735_5         # hbar^2 / (2 * amu) in MeV*fm^2

# ---- reaction (d + 16C) -------------------------------------------------
BASE = dict(m_p=2.01410, z_p=1,          # deuteron
            m_t=16.01470, z_t=6,         # 16C  (FRESCO masses)
            elab=23.55)                  # MeV, deuteron lab energy

# Geometry dict = ABSOLUTE Woods-Saxon depths(MeV)/radii(fm)/diffuseness(fm):
#   V,R_V,a_V  real volume;  W,R_W,a_W  imag volume;  WD,R_D,a_D  imag surface;
#   R_C  Coulomb (uniform sphere).  All radii absolute, not reduced.

def adwa_geom(p):
    """Johnson-Soper ADWA d+16C (FRESCO kp=2); radii = r * A_t^(1/3)."""
    a13 = p["m_t"] ** (1.0 / 3.0)
    return dict(V=102.108, R_V=1.1430 * a13, a_V=0.6754,
                W=2.030,   R_W=1.1430 * a13, a_W=0.6754,
                WD=15.540, R_D=1.3025 * a13, a_D=0.5345,
                R_C=1.4357 * a13)


def da1p_geom(p):
    """DA1p global deuteron OMP for 1p-shell / light targets
    (Pang et al., arXiv:1606.01507; Table II 1p-shell column).  Imag volume
    and surface share the geometry (R_W, a_W)."""
    A = p["m_t"]; A13 = A ** (1.0 / 3.0); E = p["elab"]
    zz_e2 = p["z_p"] * p["z_t"] * E2
    Vr, Ve = 98.9, -0.278
    rr, rr0, rre, ar = 1.11, -0.172, 0.00117, 0.776
    Wv0, Ws0 = 11.5, 7.56
    rw, rw0, rwe, aw = 0.561, 3.07, -0.00449, 0.744
    Wve0, Wvew, Wse0, Wsew = 18.1, 5.97, 14.3, 4.55
    rc = 1.3
    R_C = rc * A13
    EC = 6.0 * zz_e2 / (5.0 * R_C)          # Coulomb correction to E
    Em = E - EC
    V = Vr + Ve * Em
    R_V = rr * A13 + rr0 + rre * Em
    W = Wv0 / (1.0 + math.exp((Wve0 - Em) / Wvew))
    WD = Ws0 / (1.0 + math.exp((Em - Wse0) / Wsew))
    R_W = rw * A13 + rw0 + rwe * Em
    return dict(V=V, R_V=R_V, a_V=ar, W=W, R_W=R_W, a_W=aw,
                WD=WD, R_D=R_W, a_D=aw, R_C=R_C)


def make_potential(geom, p):
    """Return U(r) [complex MeV] and R_C from an absolute-geometry dict."""
    V, R_V, a_V = geom["V"], geom["R_V"], geom["a_V"]
    W, R_W, a_W = geom["W"], geom["R_W"], geom["a_W"]
    WD, R_D, a_D = geom["WD"], geom["R_D"], geom["a_D"]
    R_C = geom["R_C"]
    zz_e2 = p["z_p"] * p["z_t"] * E2

    def U(r):
        fV = 1.0 / (1.0 + math.exp((r - R_V) / a_V))
        fW = 1.0 / (1.0 + math.exp((r - R_W) / a_W))
        xD = (r - R_D) / a_D
        eD = math.exp(xD)
        gD = 4.0 * eD / (1.0 + eD) ** 2          # surface form, peaks at WD
        if r > R_C:
            VC = zz_e2 / r
        else:
            VC = zz_e2 / (2.0 * R_C) * (3.0 - (r / R_C) ** 2)
        return VC - V * fV - 1j * (W * fW + WD * gD)

    return U, R_C


def kinematics(p):
    mu = p["m_p"] * p["m_t"] / (p["m_p"] + p["m_t"])      # reduced mass [u]
    Ecm = p["elab"] * p["m_t"] / (p["m_p"] + p["m_t"])    # non-rel CM energy
    hbar2_2mu = HBAR2_2U / mu                             # MeV*fm^2
    k = math.sqrt(Ecm / hbar2_2mu)                        # fm^-1
    eta = p["z_p"] * p["z_t"] * E2 / (2.0 * hbar2_2mu * k)
    return dict(mu=mu, Ecm=Ecm, hbar2_2mu=hbar2_2mu, k=k, eta=eta)


_COULOMB_CACHE = {}


def _coulomb_cached(L, eta, rho):
    key = (L, round(eta, 10), round(rho, 8))
    v = _COULOMB_CACHE.get(key)
    if v is None:
        v = (float(mp.coulombf(L, eta, rho)), float(mp.coulombg(L, eta, rho)),
             float(mp.coulombf(L + 1, eta, rho)), float(mp.coulombg(L + 1, eta, rho)))
        _COULOMB_CACHE[key] = v
    return v


def s_matrix(L, U, kin, h=0.02, r_match=25.0):
    """Integrate the regular radial solution and return the nuclear S_L."""
    k = kin["k"]
    inv_h2_2mu = 1.0 / kin["hbar2_2mu"]
    LL = L * (L + 1)
    N = int(round(r_match / h)) + 2

    def Wfun(r):
        return LL / (r * r) + U(r) * inv_h2_2mu - k * k     # complex

    # Numerov for u'' = W u, with f_n = 1 - (h^2/12) W_n:
    #   f_{n+1} u_{n+1} = (12 - 10 f_n) u_n - f_{n-1} u_{n-1}
    c = h * h / 12.0
    # Seed the first two nonzero points with the regular threshold behavior
    # u(r) ~ r^(L+1) (overall scale arbitrary; the log-derivative is invariant).
    samples = {1: (1.0 + 0j), 2: (2.0 ** (L + 1)) + 0j}    # r=h, r=2h
    f_m1 = 1.0 - c * Wfun(h)
    f_0 = 1.0 - c * Wfun(2.0 * h)
    u_m1, u_0 = samples[1], samples[2]
    for n in range(3, N + 1):
        r_next = n * h
        f_next = 1.0 - c * Wfun(r_next)
        u_next = ((12.0 - 10.0 * f_0) * u_0 - f_m1 * u_m1) / f_next
        # renormalize to avoid overflow in the classically-forbidden region
        if abs(u_next) > 1e120:
            s = 1e-120
            u_next *= s; u_0 *= s; u_m1 *= s
            for kk in samples:
                samples[kk] *= s
        u_m1, u_0 = u_0, u_next
        f_m1, f_0 = f_0, f_next
        samples[n] = u_next
    # log-derivative at a = r_match via 5-point stencil (O(h^4) accurate)
    Na = int(round(r_match / h))
    a = Na * h
    u_a = samples[Na]
    uprime = (samples[Na - 2] - 8.0 * samples[Na - 1]
              + 8.0 * samples[Na + 1] - samples[Na + 2]) / (12.0 * h)
    beta = uprime / u_a                              # u'/u  [fm^-1]

    rho = k * a
    eta = kin["eta"]
    # Coulomb functions depend only on (L, eta, rho) -- fixed by kinematics,
    # not by the potential -- so cache them (huge speedup when fitting).
    F_L, G_L, F_L1, G_L1 = _coulomb_cached(L, eta, rho)
    # Coulomb derivative recurrence (verified vs high-precision finite diff):
    #   X_L'(rho) = S_{L+1} X_L - R_{L+1} X_{L+1}
    #   S_{L+1} = (L+1)/rho + eta/(L+1),  R_{L+1} = sqrt(1 + (eta/(L+1))^2)
    Lp1 = L + 1.0
    Sfac = Lp1 / rho + eta / Lp1
    Rfac = math.sqrt(1.0 + (eta / Lp1) ** 2)
    dF = Sfac * F_L - Rfac * F_L1          # dF_L/drho
    dG = Sfac * G_L - Rfac * G_L1          # dG_L/drho
    Hp = G_L + 1j * F_L                    # H^+ = G + iF  (outgoing)
    Hm = G_L - 1j * F_L                    # H^- = G - iF  (incoming)
    dHp = dG + 1j * dF
    dHm = dG - 1j * dF
    Lb = beta / k                          # dimensionless log-deriv wrt rho
    S_L = (dHm - Lb * Hm) / (dHp - Lb * Hp)
    return S_L


def coulomb_sigma(Lmax, eta):
    """Coulomb phase shifts sigma_L = arg Gamma(L+1+i eta)."""
    sig = np.zeros(Lmax + 1)
    sig[0] = float(mp.arg(mp.gamma(mp.mpf(1) + 1j * eta)))
    for L in range(1, Lmax + 1):
        sig[L] = sig[L - 1] + math.atan2(eta, L)
    return sig


def cross_section(p, geom, thetas_deg, Lmax=60, h=0.02, r_match=25.0):
    U, R_C = make_potential(geom, p)
    kin = kinematics(p)
    k, eta = kin["k"], kin["eta"]
    sigma = coulomb_sigma(Lmax, eta)

    S = np.ones(Lmax + 1, dtype=complex)
    for L in range(Lmax + 1):
        S[L] = s_matrix(L, U, kin, h=h, r_match=r_match)

    th = np.radians(thetas_deg)
    cth = np.cos(th)
    sin2 = np.sin(th / 2.0) ** 2

    # Coulomb (Rutherford) amplitude
    sig0 = sigma[0]
    fC = (-eta / (2.0 * k * sin2)
          * np.exp(-1j * eta * np.log(sin2) + 2j * sig0))

    # Nuclear amplitude
    fN = np.zeros_like(th, dtype=complex)
    for L in range(Lmax + 1):
        PL = eval_legendre(L, cth)
        fN += (2 * L + 1) * cmath.exp(2j * sigma[L]) * (S[L] - 1.0) * PL
    fN *= 1.0 / (2j * k)

    f = fC + fN
    dsdo = np.abs(f) ** 2 * 10.0          # fm^2/sr -> mb/sr
    ruth = (eta / (2.0 * k * sin2)) ** 2 * 10.0
    return dsdo, ruth, kin, S


def forward_norm(ang, dat, err, om, cut=100.0):
    """Error-weighted <data/OM> over the forward lobe (both > cut mb/sr)."""
    m = (dat > cut) & (om > cut)
    w = 1.0 / err[m] ** 2
    r = dat[m] / om[m]
    norm = np.sum(w * r) / np.sum(w)
    spread = np.sqrt(np.sum(w * (r - norm) ** 2) / np.sum(w))
    return norm, spread, ang[m].min(), ang[m].max(), int(m.sum())


def best_scale(dat, err, om):
    """Single multiplicative factor N minimizing chi^2(data, N*om), all angles.
    N is 'how much the data wants to rescale the model' = norm correction."""
    w = 1.0 / err ** 2
    return np.sum(w * dat * om) / np.sum(w * om ** 2)


def chi2_ndf(dat, err, om, npar=0):
    return np.sum(((dat - om) / err) ** 2) / (len(dat) - npar)


def fit_potential(p, ang, dat, err, start_geom, lmax):
    """Fit {V, W, WD, R_V, R_W} starting from start_geom to the full angular
    distribution (absolute, no free normalization)."""
    from scipy.optimize import least_squares
    g0 = dict(start_geom)
    x0 = np.array([g0["V"], g0["W"], g0["WD"], g0["R_V"], g0["R_W"]])

    def geom_of(x):
        g = dict(g0)
        g["V"], g["W"], g["WD"], g["R_V"], g["R_W"] = x
        g["R_D"] = g["R_W"]          # imag surface shares imag-volume geometry
        return g

    def resid(x):
        if np.any(x[:3] < 0) or np.any(x[3:] < 1.0):
            return np.full(len(dat), 1e3)
        om, _, _, _ = cross_section(p, geom_of(x), ang, Lmax=lmax)
        return (dat - om) / err

    sol = least_squares(resid, x0, method="lm", max_nfev=200)
    return geom_of(sol.x), sol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/Users/quantumlab/Downloads/16C_dd_gs.txt")
    ap.add_argument("--out", default="/Users/quantumlab/C17_dp_fits/results/om_elastic_dd.csv")
    ap.add_argument("--lmax", type=int, default=60)
    ap.add_argument("--no-fit", action="store_true", help="skip the potential fit")
    ap.add_argument("--selftest", action="store_true",
                    help="point Coulomb, no nuclear: output must equal Rutherford")
    args = ap.parse_args()

    p = dict(BASE)
    data = np.loadtxt(args.data, skiprows=1)
    ang, dat, err = data[:, 0], data[:, 1], data[:, 2]

    if args.selftest:
        geom = adwa_geom(p)
        geom.update(V=0.0, W=0.0, WD=0.0, R_C=1e-4)
        dsdo, ruth, kin, _ = cross_section(p, geom, ang, Lmax=args.lmax)
        rel = np.abs(dsdo - ruth) / ruth
        print(f"# SELFTEST max|OM-Ruth|/Ruth = {rel.max():.2e}  (should be ~1e-3)")
        return

    kin = kinematics(p)
    print(f"# d+16C elastic OM   Elab={p['elab']} MeV  Ecm={kin['Ecm']:.3f} MeV  "
          f"k={kin['k']:.4f} fm^-1  eta={kin['eta']:.4f}")

    # --- the three models ----
    models = []   # (label, geom, color, npar)
    g_adwa = adwa_geom(p); g_da1p = da1p_geom(p)
    om_adwa, ruth, _, _ = cross_section(p, g_adwa, ang, Lmax=args.lmax)
    om_da1p, _, _, _ = cross_section(p, g_da1p, ang, Lmax=args.lmax)
    models.append(("ADWA (transfer OMP)", g_adwa, om_adwa, "crimson", 0))
    models.append(("DA1p (light-target OMP)", g_da1p, om_da1p, "royalblue", 0))

    if not args.no_fit:
        print("# fitting potential (start = DA1p) ...")
        g_fit, sol = fit_potential(p, ang, dat, err, g_da1p, args.lmax)
        om_fit, _, _, _ = cross_section(p, g_fit, ang, Lmax=args.lmax)
        models.append(("fit to elastic", g_fit, om_fit, "green", 5))
        print(f"#   best-fit: V={g_fit['V']:.2f} R_V={g_fit['R_V']:.3f}  "
              f"W={g_fit['W']:.2f} WD={g_fit['WD']:.2f} R_W={g_fit['R_W']:.3f}")

    # --- report per model ----
    print(f"\n{'model':26s} {'chi2/ndf':>9s} {'fwd<d/OM>':>11s} "
          f"{'best-scale N':>13s}")
    summary = []
    for label, geom, om, color, npar in models:
        norm, spread, lo, hi, n = forward_norm(ang, dat, err, om)
        N = best_scale(dat, err, om)
        c2 = chi2_ndf(dat, err, om, npar)
        print(f"{label:26s} {c2:9.2f} {norm:6.3f}+/-{spread:.3f} {N:13.3f}")
        summary.append((label, c2, norm, spread, N))

    # CSV with all model curves
    cols = [ang, dat, err, ruth] + [m[2] for m in models]
    head = "ang_deg data_mb err_mb Ruth_mb " + " ".join(
        m[0].split()[0] + "_mb" for m in models)
    np.savetxt(args.out, np.column_stack(cols), header=head, fmt="%12.5f")
    print(f"\n# wrote {args.out}")

    # Interpretation: the forward Coulomb-nuclear lobe is where the absolute
    # scale is robust against the potential.  Headline = the fitted potential
    # (best chi^2 + best shape); spread across potentials = OM-model systematic.
    fwd_all = [s[2] for s in summary]
    ref_label, _, ref_norm, ref_sp, _ = summary[-1]   # fit if present, else DA1p
    print(f"\n# forward-lobe <data/OM>:  {', '.join(f'{l.split()[0]}={n:.3f}' for l,_,n,_,_ in summary)}")
    print(f"# headline ({ref_label}) = {ref_norm:.3f} +/- {ref_sp:.3f}; "
          f"cross-potential spread {min(fwd_all):.3f}-{max(fwd_all):.3f}")
    if abs(ref_norm - 1.0) < 0.15:
        print(f"# => normalization consistent with 1.0 within ~{abs(ref_norm-1)*100:.0f}%: "
              "(d,p) ABSOLUTE NORMALIZATION VALIDATED (no rescale needed).")
    else:
        print(f"# => (d,p) cross sections would need scaling by {ref_norm:.3f}.")

    # ---- overlay plot ---------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 5.7))
        ax.errorbar(ang, dat, yerr=err, fmt="o", ms=4, color="k",
                    capsize=2, label=r"d+$^{16}$C elastic (data)", zorder=5)
        ax.plot(ang, ruth, ":", lw=1.1, color="0.55", label="Rutherford")
        for label, geom, om, color, npar in models:
            ax.plot(ang, om, "-", lw=1.8, color=color, label=label)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\theta_{\rm c.m.}$ (deg)", fontsize=13)
        ax.set_ylabel(r"d$\sigma$/d$\Omega$ (mb/sr)", fontsize=13)
        ax.set_title(r"$^{16}$C(d,d)$^{16}$C elastic, $E_{\rm lab}$=23.55 MeV "
                     "-- absolute-normalization check", fontsize=11)
        ax.set_xlim(10, 91); ax.set_ylim(1e-2, 3e3)
        ax.legend(fontsize=9.5, frameon=False)
        ax.tick_params(top=True, right=True, direction="in", which="both")
        fig.tight_layout()
        plot_path = args.out.replace("results/", "plots/").replace(".csv", ".png")
        fig.savefig(plot_path, dpi=140)
        print(f"# wrote {plot_path}")
    except Exception as e:  # pragma: no cover
        print(f"# (plot skipped: {e})")


if __name__ == "__main__":
    main()
