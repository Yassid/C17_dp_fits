// C16_pd_AngBins.C  --  angle-binned spectral fitting for 16C(d,p)17C
//
// Refactor of C16_pd_ana_penetrability_ang_dist.C with:
//   * fixed local paths (data in /home/yassid/C16_dp/C16_dp/InterpSolver_root/,
//     penetrability and phase-space files relative to Codes/)
//   * 6 theta_CM bins fitted in one run (10-15, 15-20, 20-25, 25-30, 30-35, 35-40 deg)
//   * inclusive fit pins BW positions, widths and resolution; per-bin fits free
//     only the amplitudes (Gaussian + BW + PS + linear bg)
//   * ConvolutedBW replaced by analytic TMath::Voigt (Gamma_eff is constant inside
//     the convolution by construction, so the integrand is exactly a Voigt)
//   * yields per state per bin and dsigma/dOmega tables written to CSV
//
// Run from /home/yassid/C16dp_fits/Codes/ after
//   source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh
// then
//   root -l -b -q 'C16_pd_AngBins.C+'

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/MinimizerOptions.h"

#include "../Penetrabilities_def/penetrabilities_neutron_17C_def_L_0.C"
#include "../Penetrabilities_def/penetrabilities_neutron_17C_def_L_1.C"
#include "../Penetrabilities_def/penetrabilities_neutron_17C_def_L_2.C"
#include "../Penetrabilities_def/penetrabilities_neutron_17C_def_L_3.C"

// ============================== constants ===================================

namespace {
constexpr double kEbeam = 188.33;            // total MeV
constexpr double kU2MeV = 931.49401;
constexpr double kM_p   = 1.007825      * kU2MeV;
constexpr double kM_d   = 2.0135532     * kU2MeV;
constexpr double kM_C16 = 16.0147       * kU2MeV;
constexpr double kM_C17 = 17.0226       * kU2MeV;
constexpr int    kZ_ej  = 1;

// Excitation-axis binning -- 0.1 MeV/bin to match the github hexCorr
// (600 bins in [-5,55]) and the PDF presentation density.
constexpr double kEbinMin = -1.5;
constexpr double kEbinMax =  7.0;       // fit upper limit (1n threshold is ~0.735)
constexpr int    kNumberBins = 85;      // (kEbinMax - kEbinMin) / 0.1 = 85
constexpr double kQcorrShift = 0.208;   // physical shift applied after the per-event QcorrZ

// theta_CM bin grid: 12 bins of 2.5 deg from 10 to 40 deg, matching the
// hex11..hex43 layout of C16_pd_ana_penetrability_ang_dist.C in the github
// repo (storage binning).  Efficiency is interpolated linearly between the
// 5-deg table entries.
struct AngBin { double lo; double hi; double center() const { return 0.5 * (lo + hi); } };
const std::array<AngBin, 12> kBins = {{
    {10.0,12.5},{12.5,15.0},{15.0,17.5},{17.5,20.0},
    {20.0,22.5},{22.5,25.0},{25.0,27.5},{27.5,30.0},
    {30.0,32.5},{32.5,35.0},{35.0,37.5},{37.5,40.0}
}};

const std::array<const char*, 12> kBinTag = {
    "10-12.5","12.5-15","15-17.5","17.5-20",
    "20-22.5","22.5-25","25-27.5","27.5-30",
    "30-32.5","32.5-35","35-37.5","37.5-40"
};

// State labels and L assignments matching SpectralModel.
const std::array<const char*, 3> kGaussLabel = {
    "gs_3half+", "ex_1half+_217", "ex_5half+_335"
};
const std::array<const char*, 7> kBWLabel = {
    "ex2.763_1half-", "ex2.980_3half+", "ex3.661_NEW",
    "ex4.231_3half+", "ex4.841_1half+", "ex5.91_3half+", "ex6.30_5half+"
};
const std::array<int, 7> kBW_L = {1, 2, 2, 2, 2, 2, 2};
}  // namespace

// ============================== utilities ===================================

double omega(double x, double y, double z) {
    return std::sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z);
}

std::tuple<double, double> kine_2b(double m1, double m2, double m3, double m4,
                                   double K_proj, double thetalab, double K_eject) {
    double Et1 = K_proj + m1;
    double Et2 = m2;
    double Et3 = K_eject + m3;
    double Et4 = Et1 + Et2 - Et3;
    double s = m1*m1 + m2*m2 + 2*m2*Et1;
    double u = m2*m2 + m3*m3 - 2*m2*Et3;
    double m4_ex = std::sqrt((std::cos(thetalab) * omega(s, m1*m1, m2*m2) * omega(u, m2*m2, m3*m3) -
                              (s - m1*m1 - m2*m2) * (m2*m2 + m3*m3 - u)) / (2*m2*m2) +
                             s + u - m2*m2);
    double Ex = m4_ex - m4;
    double t  = m2*m2 + m4_ex*m4_ex - 2*m2*Et4;
    double theta_cm = TMath::Pi() - std::acos(
        (s*s + s*(2*t - m1*m1 - m2*m2 - m3*m3 - m4_ex*m4_ex) +
         (m1*m1 - m2*m2) * (m3*m3 - m4_ex*m4_ex)) /
        (omega(s, m1*m1, m2*m2) * omega(s, m3*m3, m4_ex*m4_ex)));
    return {Ex, theta_cm * TMath::RadToDeg()};
}

TH1F* ShiftHistogramX(const TH1F* h, double shift) {
    int nbins = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin() + shift;
    double xmax = h->GetXaxis()->GetXmax() + shift;
    auto* hshift = new TH1F(TString(h->GetName()) + "_shift", h->GetTitle(),
                            nbins, xmin, xmax);
    for (int i = 1; i <= nbins; ++i) {
        hshift->SetBinContent(i, h->GetBinContent(i));
        hshift->SetBinError(i, h->GetBinError(i));
    }
    return hshift;
}

TGraph* histoToTgraph(TH1F* h, const char* name = "g_PS") {
    auto* g = new TGraph();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        g->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
    }
    g->SetName(name);
    return g;
}

// ============================== spectral model ==============================

// BW(t; E0, Gamma)  =  Gamma / [(t - E0)^2 + (Gamma/2)^2]
// Area over t equals 2*pi.  ConvolutedBW returns Amp * 2*pi * Voigt
// (unit-area Voigt centered at E_eff with widths Gamma_eff and sigma).
double ConvolutedBW(double E, double Amp, double E_eff, double Gamma_eff, double sigma) {
    return Amp * 2.0 * TMath::Pi() *
           TMath::Voigt(E - E_eff, sigma, Gamma_eff, 4);
}

// Penetrability-modified BW.  Pen[L] tables come from the included headers.
double BWModificada(double E, double Amp, double E0, double Gamma0, double sigma,
                    int l) {
    constexpr double E_min = 0.735;     // 1n threshold (n + 16C)
    constexpr double E_max = 7.0;
    constexpr int    num_bins_bw = 10000;
    constexpr int    num_bins_pen = 10000;
    constexpr double dE_bw  = 0.0006265;
    const     double dE_pen = (E_max - E_min) / num_bins_pen;

    int bin_index_bw = static_cast<int>((E - E_min) / dE_bw);
    if (bin_index_bw < 0 || bin_index_bw >= num_bins_bw) return 0.0;

    int bin_index_E0 = static_cast<int>((E0 - E_min) / dE_pen);
    bin_index_E0 = std::clamp(bin_index_E0, 0, num_bins_pen - 1);

    int bin_index_pen = static_cast<int>((E - E_min) / dE_pen);
    bin_index_pen = std::min(bin_index_pen, num_bins_pen - 1);

    double E_bin = E_min + (bin_index_bw + 0.5) * dE_bw;

    double* T[4] = {
        T0_neutron_17C_values, T1_neutron_17C_values,
        T2_neutron_17C_values, T3_neutron_17C_values
    };
    double Gamma_eff = Gamma0;
    double E_eff     = E0;
    if (E_bin >= E_min) {
        Gamma_eff = Gamma0 * T[l][bin_index_pen] / T[l][bin_index_E0];
    }
    return ConvolutedBW(E, Amp, E_eff, Gamma_eff, sigma);
}

class SpectralModel {
public:
    TGraph* g_PS1n;
    TGraph* g_PS2n;
    SpectralModel(TGraph* g, TGraph* g2) : g_PS1n(g), g_PS2n(g2) {}

    double operator()(double* x, double* p) {
        double val = 0;
        // 3 Gaussians (bound states): amp, mean, sigma triples at p[0..8]
        for (int g = 0; g < 3; ++g) {
            val += p[3*g] * TMath::Gaus(x[0], p[3*g + 1], p[3*g + 2], false);
        }
        // 7 BWs at p[9..36]: each (Amp, E0, Gamma0, sigma)
        for (int b = 0; b < 7; ++b) {
            int o = 9 + 4*b;
            val += BWModificada(x[0], p[o], p[o + 1], p[o + 2], p[o + 3], kBW_L[b]);
        }
        // p[37] = 1n PS amplitude, p[38] = linear bg, p[39] = 2n PS amplitude
        val += p[37] * g_PS1n->Eval(x[0]);
        val += p[38];
        val += p[39] * g_PS2n->Eval(x[0]);
        return val;
    }
};

// ============================== runtime =====================================

void C16_pd_AngBins() {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    gStyle->SetOptStat(0);

    // ---------------- data loading ----------------
    const TString dataDir = "/home/yassid/C16_dp/C16_dp/InterpSolver_root/";
    // Full run list: union of the parent C16_pd_ana.C (45 runs starting at
    // 0016) and the github _ang_dist.C macro (28 runs ending at 0103).
    // All 47 are present in InterpSolver_root/.  Matching this set is what
    // recovers the higher statistics shown in the PDF presentation.
    std::vector<TString> filenames = {
        "run_0016_1H.root","run_0017_1H.root","run_0018_1H.root","run_0019_1H.root",
        "run_0020_1H.root","run_0021_1H.root","run_0022_1H.root","run_0023_1H.root",
        "run_0026_1H.root","run_0027_1H.root","run_0031_1H.root","run_0032_1H.root",
        "run_0034_1H.root","run_0036_1H.root","run_0037_1H.root","run_0038_1H.root",
        "run_0039_1H.root","run_0040_1H.root","run_0041_1H.root","run_0042_1H.root",
        "run_0043_1H.root","run_0044_1H.root","run_0046_1H.root","run_0048_1H.root",
        "run_0057_1H.root","run_0058_1H.root","run_0076_1H.root","run_0077_1H.root",
        "run_0078_1H.root","run_0079_1H.root","run_0080_1H.root","run_0082_1H.root",
        "run_0083_1H.root","run_0084_1H.root","run_0085_1H.root","run_0086_1H.root",
        "run_0087_1H.root","run_0088_1H.root","run_0089_1H.root","run_0091_1H.root",
        "run_0092_1H.root","run_0095_1H.root","run_0096_1H.root","run_0097_1H.root",
        "run_0098_1H.root","run_0102_1H.root","run_0103_1H.root"
    };

    // Inclusive (10-40) and per-bin Ex histograms.  We fill Q-corrected Ex.
    auto* hInclusive = new TH1F("hInclusive", "Inclusive 10-40 deg",
                                kNumberBins, kEbinMin, kEbinMax);
    std::array<TH1F*, 12> hBin;
    for (size_t b = 0; b < kBins.size(); ++b) {
        hBin[b] = new TH1F(Form("hBin_%zu", b),
                           Form("17C Ex spectrum theta_CM %s deg", kBinTag[b]),
                           kNumberBins, kEbinMin, kEbinMax);
    }
    auto* hThetaCM = new TH1F("hThetaCM", "theta_CM all events", 90, 0, 90);

    // Per-event Q-correction constants kept identical to the existing analysis.
    constexpr double kP0 = 0.618451;
    constexpr double kP1 = 0.0080509149;
    constexpr double kOffset = 1.022;

    long total_events = 0, kept_events = 0;
    for (auto& filename : filenames) {
        TFile* runFile = TFile::Open(dataDir + filename, "READ");
        if (!runFile || runFile->IsZombie()) {
            std::cerr << "  skip (missing): " << filename << "\n";
            continue;
        }
        auto* Tphysics = static_cast<TTree*>(runFile->Get("parquettree"));
        double theta = 0, phi = 0, Brho = 0, redchi = 0, zPos = 0;
        Tphysics->SetBranchAddress("polar", &theta);
        Tphysics->SetBranchAddress("azimuthal", &phi);
        Tphysics->SetBranchAddress("brho", &Brho);
        Tphysics->SetBranchAddress("redchisq", &redchi);
        Tphysics->SetBranchAddress("vertex_z", &zPos);

        const Long64_t n_ev = Tphysics->GetEntries();
        for (Long64_t i = 0; i < n_ev; ++i) {
            Tphysics->GetEntry(i);
            ++total_events;
            const double p_ej = Brho * kZ_ej * 2.99792458e2;          // MeV/c
            const double E_ej = std::sqrt(p_ej*p_ej + kM_p*kM_p) - kM_p;
            // Cut set matching the parent C16_pd_ana.C (used for the PDF
            // presentation): only the polar-angle floor.  The github
            // _ang_dist.C added E_ej<8 + zPos in [2,98] which discard ~3x
            // events and were not used for the PDF figures.
            const double thetaDeg = theta * TMath::RadToDeg();
            if (thetaDeg < 90.0) continue;
            const double zCm = zPos * 100.0;

            auto [ex_raw, theta_cm] = kine_2b(kM_C16, kM_d, kM_p, kM_C17,
                                              kEbeam, theta, E_ej);
            const double Qcorr = ex_raw - kP1 * zCm - kP0 + kOffset;
            const double Ex_corr = Qcorr + kQcorrShift;

            hThetaCM->Fill(theta_cm);

            // Inclusive 10-40
            if (theta_cm >= kBins.front().lo && theta_cm < kBins.back().hi) {
                hInclusive->Fill(Ex_corr);
                ++kept_events;
                for (size_t b = 0; b < kBins.size(); ++b) {
                    if (theta_cm >= kBins[b].lo && theta_cm < kBins[b].hi) {
                        hBin[b]->Fill(Ex_corr);
                        break;
                    }
                }
            }
        }
        runFile->Close();
        delete runFile;
    }
    std::cout << "Events: total " << total_events << "  kept in 10-40 deg " << kept_events << "\n";
    for (size_t b = 0; b < kBins.size(); ++b) {
        std::cout << "  bin " << kBinTag[b] << ": "
                  << static_cast<long>(hBin[b]->GetEntries()) << " events\n";
    }

    // ---------------- phase-space templates (per bin + inclusive) -----------
    auto buildPS = [&](const TString& fname, double th_lo, double th_hi,
                       const char* hname) -> TH1F* {
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Phase-space file not found: " << fname << "\n";
            return new TH1F(hname, hname, kNumberBins, kEbinMin, kEbinMax);
        }
        auto* t = static_cast<TTree*>(f->Get("simulated_tree"));
        double Weight_sim = 0, Ex_cal = 0, ThetaCM_cal = 0;
        t->SetBranchAddress("Weight_sim", &Weight_sim);
        t->SetBranchAddress("Ex_cal", &Ex_cal);
        t->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal);
        auto* h = new TH1F(hname, hname, kNumberBins, kEbinMin, kEbinMax);
        h->SetDirectory(nullptr);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (ThetaCM_cal >= th_lo && ThetaCM_cal < th_hi) {
                h->Fill(Ex_cal, Weight_sim);
            }
        }
        h->Smooth();
        f->Close();
        return h;
    };

    const TString ps1nFile = "../Phase_Space/PhaseSpace_16C_dp_1n.root";
    const TString ps2nFile = "../Phase_Space/PhaseSpace_16C_dp_2n.root";

    TH1F* h_PS_1n_incl = buildPS(ps1nFile, kBins.front().lo, kBins.back().hi, "h_PS_1n_incl");
    TH1F* h_PS_2n_incl = buildPS(ps2nFile, kBins.front().lo, kBins.back().hi, "h_PS_2n_incl");
    std::array<TH1F*, 12> h_PS_1n_bin, h_PS_2n_bin;
    for (size_t b = 0; b < kBins.size(); ++b) {
        h_PS_1n_bin[b] = buildPS(ps1nFile, kBins[b].lo, kBins[b].hi,
                                 Form("h_PS_1n_bin%zu", b));
        h_PS_2n_bin[b] = buildPS(ps2nFile, kBins[b].lo, kBins[b].hi,
                                 Form("h_PS_2n_bin%zu", b));
    }

    // ---------------- spectral models -----------------
    auto* g_PS1n_incl = histoToTgraph(h_PS_1n_incl, "g_PS1n_incl");
    auto* g_PS2n_incl = histoToTgraph(h_PS_2n_incl, "g_PS2n_incl");
    SpectralModel* mInclusive = new SpectralModel(g_PS1n_incl, g_PS2n_incl);
    TF1* fInclusive = new TF1("fInclusive", mInclusive, kEbinMin, kEbinMax, 40, "SpectralModel");

    // Initial values: amplitudes scaled by entries, parameters from the PDF table.
    const double scaleIncl = std::max(1.0, hInclusive->GetMaximum());
    std::vector<double> p0(40, 0.0);
    // Gaussians (gs, 0.218, 0.350): (Amp, mean, sigma)
    p0[0] = 0.4 * scaleIncl;  p0[1] = 0.000; p0[2] = 0.120;
    p0[3] = 0.6 * scaleIncl;  p0[4] = 0.218; p0[5] = 0.120;
    p0[6] = 0.2 * scaleIncl;  p0[7] = 0.350; p0[8] = 0.120;
    // BW: (Amp, E0, Gamma0, sigma) -- start from PDF centroids and widths
    const double bwE0[7]   = {2.763, 2.980, 3.661, 4.231, 4.841, 5.910, 6.300};
    const double bwGam0[7] = {0.050, 0.200, 0.398, 0.500, 0.300, 1.500, 0.396};
    const double bwAmp0[7] = {  20,   40,    80,    80,    40,   100,    80};
    for (int b = 0; b < 7; ++b) {
        int o = 9 + 4*b;
        p0[o]   = bwAmp0[b];
        p0[o+1] = bwE0[b];
        p0[o+2] = bwGam0[b];
        p0[o+3] = 0.120;
    }
    p0[37] = 0.010;   // 1n PS amplitude
    p0[38] = 1.000;   // linear bg constant
    p0[39] = 0.500;   // 2n PS amplitude
    for (size_t i = 0; i < p0.size(); ++i) fInclusive->SetParameter(i, p0[i]);

    // Gaussian limits/locks
    fInclusive->SetParLimits(0, 1, 5*scaleIncl);   fInclusive->FixParameter(1, 0.000); fInclusive->FixParameter(2, 0.120);
    fInclusive->SetParLimits(3, 1, 5*scaleIncl);   fInclusive->FixParameter(4, 0.218); fInclusive->FixParameter(5, 0.120);
    fInclusive->SetParLimits(6, 1, 5*scaleIncl);   fInclusive->FixParameter(7, 0.350); fInclusive->FixParameter(8, 0.120);
    // BW limits (centroid +/- a few sigma, Gamma in PDF range)
    const double bwCenLo[7] = {2.65, 2.98, 3.50, 4.10, 4.80, 5.85, 6.30};
    const double bwCenHi[7] = {2.80, 3.02, 3.70, 4.40, 5.00, 6.10, 6.50};
    const double bwGamLo[7] = {0.025, 0.10, 0.10, 0.10, 0.10, 0.25, 0.15};
    const double bwGamHi[7] = {0.080, 0.30, 0.50, 0.70, 0.50, 1.80, 0.55};
    for (int b = 0; b < 7; ++b) {
        int o = 9 + 4*b;
        fInclusive->SetParLimits(o,     0.1,  10*scaleIncl);
        fInclusive->SetParLimits(o + 1, bwCenLo[b], bwCenHi[b]);
        fInclusive->SetParLimits(o + 2, bwGamLo[b], bwGamHi[b]);
        fInclusive->FixParameter(o + 3, 0.120);
    }
    fInclusive->SetParLimits(37, 0.001, 0.10);
    fInclusive->SetParLimits(38, 0.25,  3.0);
    fInclusive->SetParLimits(39, 0.05,  2.0);

    std::cout << "\n--- inclusive fit ---\n";
    hInclusive->Fit(fInclusive, "RQM0");
    hInclusive->Fit(fInclusive, "RM");      // second pass with full Minuit output

    // Snapshot positions and widths from inclusive (these will be frozen per bin)
    std::vector<double> incl(40);
    for (int i = 0; i < 40; ++i) incl[i] = fInclusive->GetParameter(i);

    // ---------------- per-bin fits ------------------
    std::array<TF1*, 12> fBin;
    std::array<SpectralModel*, 12> mBin;
    std::vector<std::vector<double>> binPars(kBins.size(), std::vector<double>(40, 0.0));
    std::vector<double>              binChi2(kBins.size(), 0.0);
    std::vector<int>                 binNdf(kBins.size(), 0);

    for (size_t b = 0; b < kBins.size(); ++b) {
        auto* g1 = histoToTgraph(h_PS_1n_bin[b], Form("g_PS1n_b%zu", b));
        auto* g2 = histoToTgraph(h_PS_2n_bin[b], Form("g_PS2n_b%zu", b));
        mBin[b] = new SpectralModel(g1, g2);
        fBin[b] = new TF1(Form("fBin_%zu", b), mBin[b], kEbinMin, kEbinMax, 40, "SpectralModel");

        // Scale amplitudes by the ratio of integrals (rough warm start)
        const double scaleBin = std::max(1.0, hBin[b]->GetMaximum());
        const double scaleRatio = scaleBin / std::max(1.0, hInclusive->GetMaximum());
        for (int i = 0; i < 40; ++i) fBin[b]->SetParameter(i, incl[i]);
        for (int g = 0; g < 3; ++g)
            fBin[b]->SetParameter(3*g, std::max(0.5, incl[3*g] * scaleRatio));
        for (int bw = 0; bw < 7; ++bw)
            fBin[b]->SetParameter(9 + 4*bw, std::max(0.5, incl[9 + 4*bw] * scaleRatio));
        fBin[b]->SetParameter(37, std::max(1e-4, incl[37] * scaleRatio));
        fBin[b]->SetParameter(38, std::max(0.10, incl[38] * scaleRatio));
        fBin[b]->SetParameter(39, std::max(1e-2, incl[39] * scaleRatio));

        // Lock everything from inclusive except amplitudes + bg
        // Gaussians: amplitude free, mean/sigma fixed
        for (int g = 0; g < 3; ++g) {
            fBin[b]->SetParLimits(3*g, 0.01, 10*scaleBin);
            fBin[b]->FixParameter(3*g + 1, incl[3*g + 1]);
            fBin[b]->FixParameter(3*g + 2, incl[3*g + 2]);
        }
        // BWs: amplitude free, E0/Gamma/sigma fixed
        for (int bw = 0; bw < 7; ++bw) {
            int o = 9 + 4*bw;
            fBin[b]->SetParLimits(o, 0.01, 20*scaleBin);
            fBin[b]->FixParameter(o + 1, incl[o + 1]);
            fBin[b]->FixParameter(o + 2, incl[o + 2]);
            fBin[b]->FixParameter(o + 3, incl[o + 3]);
        }
        fBin[b]->SetParLimits(37, 1e-5, 0.20);
        fBin[b]->SetParLimits(38, 0.05, 5.00);
        fBin[b]->SetParLimits(39, 1e-3, 3.00);

        std::cout << "\n--- per-bin fit  " << kBinTag[b] << " deg ---\n";
        hBin[b]->Fit(fBin[b], "RQM0");
        hBin[b]->Fit(fBin[b], "RQM");

        binChi2[b] = fBin[b]->GetChisquare();
        binNdf[b]  = fBin[b]->GetNDF();
        for (int i = 0; i < 40; ++i) binPars[b][i] = fBin[b]->GetParameter(i);

        std::cout << "    chi2/NDF = " << binChi2[b] << " / " << binNdf[b]
                  << " = " << binChi2[b] / std::max(1, binNdf[b]) << "\n";
    }

    // ---------------- yields ------------------
    // Yield definitions:
    //   Gaussian g (g=0..2):  Y = amp_g * sigma_g * sqrt(2*pi) / binWidth
    //   BW b (b=0..6):        Y = TF1 integral of the *single* BW component over [E_min, E_max]
    //                              / binWidth, where the TF1 has only that BW active.
    // We report integrated counts (not divided by anything yet).
    auto bwSingle = [&](double Amp, double E0, double Gamma0, double sigma, int L) {
        // returns lambda for TF1 (so we can integrate the single state's shape)
        return [=](double* x, double*) {
            return BWModificada(x[0], Amp, E0, Gamma0, sigma, L);
        };
    };

    std::ofstream csv("../results/yields.csv");
    csv << "bin_lo,bin_hi,bin_center,state,Ex_MeV,counts,counts_err\n";

    const double binWidth = (kEbinMax - kEbinMin) / kNumberBins;
    const double sqrt2pi  = std::sqrt(2.0 * TMath::Pi());

    // 2D table per bin per state for downstream conversion
    std::vector<std::vector<double>> yields(kBins.size(),
                                            std::vector<double>(10, 0.0));
    std::vector<std::vector<double>> yieldsErr(kBins.size(),
                                               std::vector<double>(10, 0.0));

    for (size_t b = 0; b < kBins.size(); ++b) {
        // Gaussians
        for (int g = 0; g < 3; ++g) {
            const double amp  = binPars[b][3*g];
            const double sig  = binPars[b][3*g + 2];
            const double mean = binPars[b][3*g + 1];
            const double Y    = amp * sig * sqrt2pi / binWidth;
            const double dAmp = fBin[b]->GetParError(3*g);
            const double dY   = Y > 0 && amp > 0 ? Y * dAmp / amp : 0.0;
            yields[b][g]    = Y;
            yieldsErr[b][g] = dY;
            csv << kBins[b].lo << "," << kBins[b].hi << "," << kBins[b].center()
                << "," << kGaussLabel[g] << "," << mean << ","
                << Y << "," << dY << "\n";
        }
        // BWs
        for (int bw = 0; bw < 7; ++bw) {
            int o = 9 + 4*bw;
            TF1 fSingle(Form("fS_b%zu_bw%d", b, bw),
                        bwSingle(binPars[b][o], binPars[b][o + 1],
                                 binPars[b][o + 2], binPars[b][o + 3], kBW_L[bw]),
                        std::max(kEbinMin, 0.735), kEbinMax, 0);
            fSingle.SetNpx(2000);
            const double Y = fSingle.Integral(std::max(kEbinMin, 0.735), kEbinMax,
                                              1e-6) / binWidth;
            const double amp  = binPars[b][o];
            const double dAmp = fBin[b]->GetParError(o);
            const double dY   = Y > 0 && amp > 0 ? Y * dAmp / amp : 0.0;
            yields[b][3 + bw]    = Y;
            yieldsErr[b][3 + bw] = dY;
            csv << kBins[b].lo << "," << kBins[b].hi << "," << kBins[b].center()
                << "," << kBWLabel[bw] << "," << binPars[b][o + 1] << ","
                << Y << "," << dY << "\n";
        }
    }
    csv.close();
    std::cout << "\nwrote yields.csv\n";

    // ---------------- dsigma/dOmega ------------------
    // Same convention as the cell-7 absolute-norm calc used in the bound-state work:
    //   sigma(mb/sr) = counts / (Nbeam * Ntarget * dOmega * eff) * conv
    // where the calibrated product (1/(Nbeam*Ntarget))*10*1H_factor was distilled
    // into a single constant.  Use the value that reproduces the cell-7 plateau.
    // From overlay_both.py:  sig_raw = (cnt/(pi*pt))/area * 10.0/0.1
    //   pi = 161460.6557954168 (beam particles)
    //   pt = 0.019632068643898506 (target factor)
    //   the extra /0.1 was an Ex bin-width undo -- our yields are already
    //   integrated counts, so we drop the /0.1.
    constexpr double kNbeam   = 161460.6557954168;
    constexpr double kNtarget = 0.019632068643898506;
    constexpr double kMbConv  = 10.0;       // unit conversion to mb

    // Load efficiency table -- table entries are at 5-deg-spaced centers
    // (2.5, 7.5, 12.5, ...).  Linearly interpolate for the 2.5-deg bin centers.
    std::vector<double> effAng, effValVec, effErrVec;
    std::ifstream eff("/home/yassid/fair_install/16C_dp/RCS/efficiency_16Cdp_0m.txt");
    {
        std::string line;
        while (std::getline(eff, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            double a, e, de;
            if (iss >> a >> e >> de) {
                effAng.push_back(a);
                effValVec.push_back(e);
                effErrVec.push_back(de);
            }
        }
    }
    auto effInterp = [&](double theta, double& val, double& err) {
        if (effAng.empty()) { val = 1.0; err = 0.0; return; }
        if (theta <= effAng.front()) { val = effValVec.front(); err = effErrVec.front(); return; }
        if (theta >= effAng.back())  { val = effValVec.back();  err = effErrVec.back();  return; }
        auto it = std::lower_bound(effAng.begin(), effAng.end(), theta);
        size_t i = std::distance(effAng.begin(), it);
        double t = (theta - effAng[i-1]) / (effAng[i] - effAng[i-1]);
        val = (1.0 - t) * effValVec[i-1] + t * effValVec[i];
        err = (1.0 - t) * effErrVec[i-1] + t * effErrVec[i];
    };

    std::ofstream dsdo("../results/dsdo.csv");
    dsdo << "bin_lo,bin_hi,bin_center,state,Ex_MeV,dsdo_mbsr,dsdo_err_mbsr,eff,eff_err\n";

    auto stateLabel = [](int idx) {
        if (idx < 3) return kGaussLabel[idx];
        return kBWLabel[idx - 3];
    };
    auto stateEx = [&](int idx, size_t b) {
        if (idx < 3) return binPars[b][3*idx + 1];
        return binPars[b][9 + 4*(idx - 3) + 1];
    };

    for (size_t b = 0; b < kBins.size(); ++b) {
        const double thetaC = kBins[b].center();
        const double th_lo  = kBins[b].lo  * TMath::DegToRad();
        const double th_hi  = kBins[b].hi  * TMath::DegToRad();
        const double dOmega = (std::cos(th_lo) - std::cos(th_hi)) * 2.0 * TMath::Pi();
        double effVal = 1.0, effErr = 0.0;
        effInterp(thetaC, effVal, effErr);
        for (int s = 0; s < 10; ++s) {
            const double Y  = yields[b][s];
            const double dY = yieldsErr[b][s];
            const double sigma_raw = (Y / (kNbeam * kNtarget)) / dOmega * kMbConv;
            const double dSigma_stat = (Y > 0)
                ? sigma_raw * std::sqrt(1.0/Y + std::pow(dY/Y, 2)) : 0.0;
            const double sigma = sigma_raw / std::max(effVal, 1e-6);
            const double dSigma = std::sqrt(std::pow(dSigma_stat / std::max(effVal, 1e-6), 2) +
                                            std::pow(sigma * effErr / std::max(effVal, 1e-6), 2));
            dsdo << kBins[b].lo << "," << kBins[b].hi << "," << thetaC << ","
                 << stateLabel(s) << "," << stateEx(s, b) << ","
                 << sigma << "," << dSigma << ","
                 << effVal << "," << effErr << "\n";
        }
    }
    dsdo.close();
    std::cout << "wrote dsdo.csv\n";

    // ---------------- plots ----------------
    // Draws data + total fit + each individual component (3 Gaussians, 7 BWs,
    // 1n PS, 2n PS, linear bg).  Colours follow the PDF presentation:
    //   data       = black markers
    //   total fit  = red, thick
    //   Gaussians  = blue
    //   BWs        = black
    //   1n PS      = green
    //   2n PS      = dark green, dashed
    //   bg         = gray, dotted
    auto drawWithComponents = [&](TH1F* h, const std::vector<double>& p,
                                  TGraph* g1, TGraph* g2,
                                  const char* tag, const char* title,
                                  bool drawLegend) -> void {
        h->SetLineColor(kBlack); h->SetMarkerStyle(20); h->SetMarkerSize(0.7);
        h->GetXaxis()->SetRangeUser(-1.0, 6.9);
        h->SetTitle(title);
        h->Draw("E");

        // Total
        auto* mTot = new SpectralModel(g1, g2);
        auto* fTot = new TF1(Form("fTot_%s", tag), mTot,
                             kEbinMin, kEbinMax, 40, "SpectralModel");
        for (int i = 0; i < 40; ++i) fTot->SetParameter(i, p[i]);
        fTot->SetLineColor(kRed); fTot->SetLineWidth(2); fTot->SetNpx(1000);
        fTot->Draw("same");

        // 3 Gaussians (blue)
        for (int g = 0; g < 3; ++g) {
            double a = p[3*g], m = p[3*g+1], s = p[3*g+2];
            auto* fG = new TF1(Form("fG_%s_%d", tag, g), "gaus",
                               kEbinMin, kEbinMax);
            fG->SetParameters(a, m, s);
            fG->SetLineColor(kBlue); fG->SetLineWidth(1); fG->SetNpx(500);
            fG->Draw("same");
        }

        // 7 BWs (black, thin)
        for (int b = 0; b < 7; ++b) {
            int o = 9 + 4*b;
            double amp = p[o], e0 = p[o+1], gam = p[o+2], sig = p[o+3];
            int    L   = kBW_L[b];
            auto* fBw = new TF1(Form("fBw_%s_%d", tag, b),
                [=](double* x, double*) {
                    return BWModificada(x[0], amp, e0, gam, sig, L);
                },
                std::max(kEbinMin, 0.735), kEbinMax, 0);
            fBw->SetLineColor(kBlack); fBw->SetLineWidth(1); fBw->SetNpx(500);
            fBw->Draw("same");
        }

        // 1n PS (green)
        double amp1 = p[37];
        auto* fPS1 = new TF1(Form("fPS1_%s", tag),
            [=](double* x, double*) { return amp1 * g1->Eval(x[0]); },
            kEbinMin, kEbinMax, 0);
        fPS1->SetLineColor(kGreen+1); fPS1->SetLineWidth(1); fPS1->SetNpx(500);
        fPS1->Draw("same");

        // 2n PS (dark green, dashed)
        double amp2 = p[39];
        auto* fPS2 = new TF1(Form("fPS2_%s", tag),
            [=](double* x, double*) { return amp2 * g2->Eval(x[0]); },
            kEbinMin, kEbinMax, 0);
        fPS2->SetLineColor(kGreen+3); fPS2->SetLineStyle(2);
        fPS2->SetLineWidth(1); fPS2->SetNpx(500);
        fPS2->Draw("same");

        // Linear bg (gray dotted)
        double bgv = p[38];
        auto* fBg = new TF1(Form("fBg_%s", tag),
            [=](double*, double*) { return bgv; },
            kEbinMin, kEbinMax, 0);
        fBg->SetLineColor(kGray+2); fBg->SetLineStyle(3); fBg->SetLineWidth(1);
        fBg->Draw("same");

        if (drawLegend) {
            auto* leg = new TLegend(0.13, 0.55, 0.40, 0.88);
            leg->SetTextSize(0.025);
            leg->SetFillStyle(0); leg->SetBorderSize(0);
            leg->AddEntry(h,    "Experimental Data",        "ep");
            leg->AddEntry(fTot, "Total Fit",                "l");
            // We need representative TF1s for the legend entries; reuse
            // a thin proxy by referencing the last-created of each.
            auto* gProxy = (TF1*)gROOT->FindObject(Form("fG_%s_0",  tag));
            auto* bwProxy= (TF1*)gROOT->FindObject(Form("fBw_%s_0", tag));
            if (gProxy)  leg->AddEntry(gProxy,  "Gaussian individual states", "l");
            if (bwProxy) leg->AddEntry(bwProxy, "BW individual states",       "l");
            leg->AddEntry(fPS1, "1n Phase Space",  "l");
            leg->AddEntry(fPS2, "2n Phase Space",  "l");
            leg->AddEntry(fBg,  "Linear background", "l");
            leg->Draw();
        }
    };

    auto* cIncl = new TCanvas("cIncl", "Inclusive 10-40 deg", 1100, 700);
    drawWithComponents(hInclusive, incl, g_PS1n_incl, g_PS2n_incl,
                       "incl", "Inclusive 10-40 deg; Excitation Energy (MeV); Counts",
                       /*drawLegend=*/true);
    cIncl->SaveAs("../plots/plots_inclusive.png");

    auto* cBin = new TCanvas("cBin", "Per-bin fits", 1800, 1300);
    cBin->Divide(4, 3);
    for (size_t b = 0; b < kBins.size(); ++b) {
        cBin->cd(b + 1);
        auto* g1 = histoToTgraph(h_PS_1n_bin[b], Form("g_PS1n_p%zu", b));
        auto* g2 = histoToTgraph(h_PS_2n_bin[b], Form("g_PS2n_p%zu", b));
        drawWithComponents(
            hBin[b], binPars[b], g1, g2, Form("b%zu", b),
            Form("17C, theta_CM %s deg; Excitation Energy (MeV); Counts", kBinTag[b]),
            /*drawLegend=*/(b == 0));
    }
    cBin->SaveAs("../plots/plots_bins.png");

    std::cout << "\nDone.\n";
}
