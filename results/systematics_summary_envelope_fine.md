# Spectral-fit systematic envelope -- 17C(d,p) unbound states

Binning: **fine**.  Variants enumerated: 82.

## Inclusive-fit chi^2/NDF grid

| cuts \ contam | nominal | + contam Gauss |
|---|---|---|
| **loose** | 1.78 | 0.84 |
| **strict** | 1.88 | 0.97 |

The contaminant Gauss at ~0.45 MeV survives both cut sets, dropping chi^2/NDF from ~1.8 (rejected) to ~0.9 (acceptable).

## Per-state yield envelope

Yield = integrated counts over theta_CM in [10, 40] deg.  Envelope spans min..max across all variants.

| state | nominal | -% | +% | comment |
|---|---|---|---|---|
| g.s. 3/2+ | **90** | **-22.6%** | **+0.2%** | bound state, contam absorbs its high-side tail |
| 2.763 1/2- | **133** | **-55.7%** | **+32.2%** | contam + BW Γ scale |
| 2.980 3/2+ | **87** | **-27.5%** | **+178.7%** | fit instability between contam and 2.980 BW |
| 3.661 NEW | **371** | **-14.6%** | **+17.3%** | modest, sensitive to contam |
| 4.231 3/2+ | **253** | **-13.9%** | **+48.3%** | BW Γ scale dominates (narrow BWs grow 4.231) |
| 4.841 1/2+ | **904** | **-27.0%** | **+10.9%** | BW Γ scale dominates (1.2 MeV state) |
| 5.91 3/2+ | **604** | **-7.5%** | **+22.5%** | modest |
| 6.30 5/2+ | 108 | -8.8% | +2.9% | robust (< 10 %) |
| contaminant (0.45 MeV) | 0 | +0.0% | +0.0% |  |

**Bold** = envelope width > 10%.

## Takeaway for the paper / slides

Reading the envelope min..max as the systematic on integrated yield, the picture splits the 7 unbound states into three groups.

**Large systematic (degeneracy with the contaminant):**
- **2.980 (3/2+)**: -28% / +179%.  The contam Gauss at 0.45 MeV overlaps with the low-energy tail of the 2.980 BW; the fit cannot fully separate them with the current parametrization. Worth flagging as a model-dependence issue, not just a statistical uncertainty.
- **2.763 (1/2-)**: -56% / +32%.  Same family of degeneracies driven by contam + BW Γ scale interplay.

**Moderate systematic (BW Γ scale dominates):**
- **4.231 (3/2+)**: -14% / +48% -- BW Γ = 0.70 (narrow) inflates it.
- **4.841 (1/2+)**: -27% / +11% -- the 1.2 MeV state is the most sensitive to the BW Γ scaling because it has the widest range.
- **5.91 (3/2+)**: -7% / +22% -- BW Γ scale.

**Robust (statistics dominate):**
- **3.661 NEW**: -15% / +17% -- modest spread despite being a new state assignment.
- **6.30 (5/2+)**: -9% / +3% -- robust across all 82 variants.
- **gs 3/2+ (bound)**: -23% / 0% -- contam steals from its high-side tail.

**Cut variation does NOT remove the contaminant.**  Strict cuts (E_ej < 8 MeV + vertex_z in [2, 98] cm) reduce statistics by ~3x but preserve the 0.45 MeV excess at essentially the same centroid; the loose cut set is therefore preferred for the final analysis and the contaminant variant is the systematic.

**Excluded variants:** the no_3661 family (a state-existence test, not a yield systematic) and the orphan contam+relaxsigma combination (a known degeneracy between two knobs covering the same excess) are removed from the envelope.
