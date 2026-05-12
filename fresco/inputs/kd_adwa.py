#!/usr/bin/env python3
"""
Compute Koning-Delaroche (2003) nucleon optical-model parameters for p+16C and
n+16C at E = E_d/2, then sum them (Johnson-Soper ADWA) to produce an effective
d+16C optical potential for FRESCO.

Reference: A.J. Koning, J.P. Delaroche, Nucl. Phys. A 713 (2003) 231.
ADWA reference: R.C. Johnson, P.J.R. Soper, Phys. Rev. C 1 (1970) 976.
"""
from math import exp


def kd_params(A: int, Z: int, E: float, projectile: str):
    """Koning-Delaroche OMP for nucleon `projectile` ('p' or 'n') on target (A,Z)
    at lab energy E (MeV). Returns dict of depths and Woods-Saxon geometry."""
    N = A - Z
    asym = (N - Z) / A
    A13 = A ** (1.0 / 3.0)
    A_m13 = A ** (-1.0 / 3.0)
    A_m23 = A ** (-2.0 / 3.0)
    A_m53 = A ** (-5.0 / 3.0)

    if projectile == "p":
        Ef = -8.4075 + 0.01378 * A
        v1 = 59.30 + 21.0 * asym - 0.024 * A
        v2 = 0.007067 + 4.23e-6 * A
        v3 = 1.729e-5 + 1.136e-8 * A
        w1 = 14.667 + 0.009629 * A
        d1 = 16.0 + 16.0 * asym
        aD = 0.5187 + 5.205e-4 * A
    elif projectile == "n":
        Ef = -11.2814 + 0.02646 * A
        v1 = 59.30 - 21.0 * asym - 0.024 * A
        v2 = 0.007228 - 1.48e-6 * A
        v3 = 1.994e-5 - 2.0e-8 * A
        w1 = 12.195 + 0.0167 * A
        d1 = 16.0 - 16.0 * asym
        aD = 0.5446 - 1.656e-4 * A
    else:
        raise ValueError("projectile must be 'p' or 'n'")

    v4 = 7e-9
    w2 = 73.55 + 0.0795 * A
    d2 = 0.0180 + 0.003802 / (1.0 + exp((A - 156.0) / 8.0))
    d3 = 11.5
    vso1 = 5.922 + 0.0030 * A
    vso2 = 0.0040
    wso1 = -3.1
    wso2 = 160.0

    rV = 1.3039 - 0.4054 * A_m13
    aV = 0.6778 - 1.487e-4 * A
    rD = 1.3424 - 0.01585 * A13
    rSO = 1.1854 - 0.647 * A_m13
    aSO = 0.59
    rC = 1.198 + 0.697 * A_m23 + 12.994 * A_m53

    EE = E - Ef
    V_V = v1 * (1.0 - v2 * EE + v3 * EE * EE - v4 * EE ** 3)
    if projectile == "p":
        VC = 1.73 * Z / (rC * A13)
        dVC = VC * v1 * (v2 - 2 * v3 * EE + 3 * v4 * EE * EE)
        V_V += dVC
    W_V = w1 * EE * EE / (EE * EE + w2 * w2)
    W_D = d1 * EE * EE * exp(-d2 * EE) / (EE * EE + d3 * d3)
    V_SO = vso1 * exp(-vso2 * EE)
    W_SO = wso1 * EE * EE / (EE * EE + wso2 * wso2)

    return dict(
        V_V=V_V, rV=rV, aV=aV,
        W_V=W_V,
        W_D=W_D, rD=rD, aD=aD,
        V_SO=V_SO, W_SO=W_SO, rSO=rSO, aSO=aSO,
        rC=rC, Ef=Ef,
    )


def fmt_pot(label, depth, r, a):
    return f"  {label}: V={depth:7.3f} MeV  r={r:.4f} fm  a={a:.4f} fm"


def main():
    A, Z = 16, 6
    Ed = 23.55
    En = Ep = Ed / 2.0

    p = kd_params(A, Z, Ep, "p")
    n = kd_params(A, Z, En, "n")

    print(f"KD OMP for target A={A}, Z={Z}, projectile energy = E_d/2 = {Ep:.3f} MeV\n")
    for label, par in [("p+16C", p), ("n+16C", n)]:
        print(f"{label}  (Ef={par['Ef']:.3f} MeV, E-Ef={Ep - par['Ef']:.3f} MeV)")
        print(fmt_pot("real volume    ", par["V_V"], par["rV"], par["aV"]))
        print(fmt_pot("imag volume    ", par["W_V"], par["rV"], par["aV"]))
        print(fmt_pot("imag surface   ", par["W_D"], par["rD"], par["aD"]))
        print(fmt_pot("real spin-orbit", par["V_SO"], par["rSO"], par["aSO"]))
        print(fmt_pot("imag spin-orbit", par["W_SO"], par["rSO"], par["aSO"]))
        print(f"  Coulomb radius rC = {par['rC']:.4f} fm")
        print()

    # Johnson-Soper sum (central terms only). Geometry: real-volume p and n
    # geometries are identical (KD depends only on A for rV, aV). For the
    # surface (W_D) we average aD slightly between p and n.
    V_d = p["V_V"] + n["V_V"]
    W_Vd = p["W_V"] + n["W_V"]
    W_Dd = p["W_D"] + n["W_D"]
    rV = p["rV"]
    aV = p["aV"]
    rD = p["rD"]
    aD = 0.5 * (p["aD"] + n["aD"])
    rC = p["rC"]

    print("Johnson-Soper ADWA effective d+16C OMP (sum of p,n central terms):")
    print(fmt_pot("real volume ", V_d, rV, aV))
    print(fmt_pot("imag volume ", W_Vd, rV, aV))
    print(fmt_pot("imag surface", W_Dd, rD, aD))
    print(f"  Coulomb radius rC (proton) = {rC:.4f} fm")
    print()

    print("FRESCO &pot lines for kp=2 (d+16C, ADWA):")
    print(f" &pot kp= 2 type= 0 p(1:3)= {A:6.3f}   0.000  {rC:6.4f} /")
    print(f" &pot kp= 2 type= 1 p(1:7)= {V_d:7.3f}  {rV:6.4f}  {aV:6.4f}  "
          f"{W_Vd:6.3f}  {rV:6.4f}  {aV:6.4f}  0.000 /")
    print(f" &pot kp= 2 type= 2 p(1:7)=  0.000   0.000   0.000  "
          f"{W_Dd:6.3f}  {rD:6.4f}  {aD:6.4f}  0.000 /")
    print(" ! spin-orbit dropped in Johnson-Soper prescription")


if __name__ == "__main__":
    main()
