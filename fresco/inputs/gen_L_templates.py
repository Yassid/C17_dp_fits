#!/usr/bin/env python3
"""
Generate four FRESCO inputs (dp17C_L0.nin, dp17C_L1.nin, dp17C_L2.nin,
dp17C_L3.nin), one per orbital angular momentum L of the transferred neutron.
Each input contains all 7 unbound resonance excitation energies as 'final
states' with the (n, l, j) configuration appropriate to that L:

  L=0 -> 3s1/2  (nn=3)
  L=1 -> 2p1/2  (nn=2)
  L=2 -> 2d3/2  (nn=2)
  L=3 -> 1f5/2  (nn=1)

All other settings (ADWA d+16C OMP, p+17C OMP, weakly-bound approx be=0.5 MeV,
FRESCO control parameters) match dp17C_unbound_adwa.nin.

After writing the files, the script optionally drives run_fresco.sh on each.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FRESCO = os.path.join(HERE, "..")  # parent of inputs/ is the fresco/ dir
RUN_SH = os.path.join(FRESCO, "run_fresco.sh")

# State Ex values (MeV) from the PDF
EX = [2.763, 2.980, 3.661, 4.231, 4.841, 5.910, 6.300]
TAGS = ["2.763", "2.980", "3.661_NEW", "4.231", "4.841", "5.91", "6.30"]

# Per-L configuration: (nn, l, j) for the transferred neutron + parity sign
# for the final-state J^π (jt, ptyt).
L_CONFIG = {
    0: dict(nn=3, l=0, j=0.5, jt=0.5, ptyt= 1, desc="3s1/2"),
    1: dict(nn=2, l=1, j=0.5, jt=0.5, ptyt=-1, desc="2p1/2"),
    2: dict(nn=2, l=2, j=1.5, jt=1.5, ptyt= 1, desc="2d3/2"),
    3: dict(nn=1, l=3, j=2.5, jt=2.5, ptyt=-1, desc="1f5/2"),
}


def make_input(L):
    cfg = L_CONFIG[L]
    nn, l, j = cfg["nn"], cfg["l"], cfg["j"]
    jt, ptyt = cfg["jt"], cfg["ptyt"]

    header = (
        f"16C(d,p)17C transfer at Elab(d)=23.55 MeV; L={L} ({cfg['desc']}) template "
        f"for 7 resonances; ADWA OMP, weakly-bound (be=0.5 MeV)\n"
        "NAMELIST\n"
        " &FRESCO  hcm= 0.05 rmatch= 60.0 rintp= 0.5\n"
        "     hnl= 0.1 rnl= 3.0 centre= 0.0 cutl=-1.6\n"
        "     jtmin= 0.0 jtmax= 40.0 absend= -1.0\n"
        "     thmin= 1.0 thmax=180.0 thinc= 2.0\n"
        "     ips= 0.05 it0=0 iter=1 iblock=0\n"
        "     chans=1 listcc=0 treneg=3 cdetr=0 smats=2 xstabl=0 veff=1\n"
        "     elab(1)= 23.550 /\n\n"
    )

    p1 = (
        " &PARTITION namep='DEUTERON' massp= 2.01410 zp= 1 nex= 1 pwf=T\n"
        "            namet='C-16    ' masst=16.01470 zt= 6 qval=-1.498/\n"
        " &STATES jp= 1.0 ptyp= 1 ep= 0.0000 cpot= 2 jt= 0.0 ptyt= 1 et= 0.0000/\n\n"
    )

    p2 = (
        " &PARTITION namep='PROTON  ' massp= 1.00728 zp= 1 nex= 7 pwf=T\n"
        "            namet='C-17    ' masst=17.02258 zt= 6 qval= 0.0000/\n"
    )
    for ex, tag in zip(EX, TAGS):
        p2 += (
            f" &STATES copyp= 1                          "
            f"jt= {jt:.1f} ptyt={ptyt:+2d} et= {ex:.4f}/  ! {tag} MeV  L={L} {cfg['desc']}\n"
        )
    # First state replaces copyp=1 with a full block to anchor jp/ep/cpot
    p2 = p2.replace(
        " &STATES copyp= 1                          ",
        " &STATES jp= 0.5 ptyp= 1 ep= 0.0000 cpot= 1 ", 1,
    )
    p2 += " &partition /\n\n"

    pot = (
        "! kp=1: p+17C optical model -- Koning-Delaroche, Ep=21 MeV, A=17, Z=6\n"
        " &pot kp= 1 type= 0 p(1:3)= 17.000   0.000   1.4190 /\n"
        " &pot kp= 1 type= 1 p(1:7)= 42.361   1.1462  0.6753  1.689   1.1462  0.6753  0.000 /\n"
        " &pot kp= 1 type= 2 p(1:7)=  0.000   0.000   0.000   5.144   1.3016  0.5275  0.000 /\n"
        " &pot kp= 1 type= 3 p(1:3)=  5.305   0.9338  0.5900 /\n\n"
        "! kp=2: d+16C Johnson-Soper ADWA from KD p+n at E_d/2 = 11.78 MeV\n"
        " &pot kp= 2 type= 0 p(1:3)= 16.000   0.000   1.4357 /\n"
        " &pot kp= 2 type= 1 p(1:7)= 102.108  1.1430  0.6754   2.030  1.1430  0.6754  0.000 /\n"
        " &pot kp= 2 type= 2 p(1:7)=  0.000   0.000   0.000   15.540  1.3025  0.5345  0.000 /\n\n"
        "! kp=3: bound-state Woods-Saxon for n+16C\n"
        " &pot kp= 3 type= 0 p(1:3)= 16.000   0.000   1.250 /\n"
        " &pot kp= 3 type= 1 p(1:3)= 50.000   1.250   0.650 /\n"
        " &pot kp= 3 type= 3 p(1:3)=  6.000   1.250   0.650 /\n\n"
        "! kp=4: deuteron internal Reid-Soft-Core potential\n"
        " &pot kp= 4 type= 0 p(1:3)=  1.000   0.000   1.250 /\n"
        " &pot kp= 4 type= 1 itt=F shape= 5 p(1:3)=  1.000   0.000   1.000 /\n"
        " &pot kp= 4 type= 3 itt=F shape= 5 p(1:3)=  1.000   0.000   1.000 /\n"
        " &pot kp= 4 type= 4 itt=F shape= 5 p(1:3)=  1.000   0.000   1.000 /\n"
        " &pot kp=-4 type= 7 itt=F shape= 5 p(1:3)=  1.000   0.000   1.000 /\n"
        " &pot /\n\n"
    )

    # 1 deuteron overlap + 7 single-particle overlaps (one per state, same L)
    overlaps = (
        " &OVERLAP kn1= 1 kn2= 2 ic1= 1 ic2= 2 in=-1 kind= 3 ch1=' ' nn= 1 l= 0 lmax= 2 sn= 0.5 ia= 0 j= 0.5 ib= 0\n"
        "    kbpot= 4 krpot= 0 be= 2.2260 isc= 0 ipc= 3 nfl= 0 nam= 0 ampl= 0.0000 /\n\n"
    )
    for i, (ex, tag) in enumerate(zip(EX, TAGS)):
        kn = 3 + i
        overlaps += (
            f"! kn={kn}: 17C({ex:.3f} {tag}) = 16C(0+) + n[{cfg['desc']}], L={l} j={j}, weakly-bound be=0.5\n"
            f" &OVERLAP kn1= {kn} kn2= 0 ic1= 1 ic2= 2 in= 2 kind= 0 nn= {nn} l= {l} sn= 0.5 j= {j:.1f}\n"
            f"    kbpot= 3 krpot= 0 be= 0.5000 isc= 1 ipc= 1 nfl= 0 nam= 0 ampl= 0.0000 /\n\n"
        )
    overlaps += " &overlap /\n\n"

    couplings = ""
    for i in range(7):
        couplings += (
            " &COUPLING icto= 2 icfrom= 1 kind= 7 ip1= 0 ip2= 0 ip3= 0 /\n"
            "   &cfp in= 1 ib= 1 ia= 1 kn= 1 a= 1.000 /\n"
            f"   &cfp in=-2 ib= {i+1} ia= 1 kn= {3+i} a= 1.000 /\n"
        )
    couplings += " &COUPLING /\n"

    return header + p1 + p2 + pot + overlaps + couplings


def main():
    files = []
    for L in (0, 1, 2, 3):
        path = os.path.join(HERE, f"dp17C_L{L}.nin")
        with open(path, "w") as f:
            f.write(make_input(L))
        print(f"wrote {path}")
        files.append(os.path.basename(path).replace(".nin", ""))

    if "--run" in sys.argv:
        for stem in files:
            print(f"\n=== running FRESCO on {stem} ===")
            r = subprocess.run([RUN_SH, stem], capture_output=True, text=True)
            tail = "\n".join(r.stdout.splitlines()[-8:])
            print(tail)
            if r.returncode != 0:
                print("STDERR:\n", r.stderr[-2000:])


if __name__ == "__main__":
    main()
