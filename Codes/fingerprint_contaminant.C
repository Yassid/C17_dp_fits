// Physical fingerprint of the 0.45 MeV "contaminant" -- self-contained
// (no #include of the analysis macro; constants + kinematics copied verbatim so
// it can run interpreted: root -l -b -q 'fingerprint_contaminant.C').
//
// Reconstructs the same Q-corrected Ex and theta_cm as C16_pd_AngBins.C for the
// 47 runs (nominal cut polar>=90), then compares the contaminant window
// (Ex in [0.30,1.10] MeV) at BACK angles (theta_cm>=30, where it lives) vs the
// same window FORWARD and vs a genuine-transfer reference (Ex in [3.5,5.0]) at
// back angles, in vertex_z / lab angle / proton energy, plus 2D Ex maps.

#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"

namespace {
constexpr double kU2MeV = 931.49410242;
constexpr double kEbeam = 188.33;
constexpr double kM_p   = 1.007825   * kU2MeV;
constexpr double kM_d   = 2.0135532  * kU2MeV;
constexpr double kM_C16 = 16.0147    * kU2MeV;
constexpr double kM_C17 = 17.0226    * kU2MeV;
constexpr int    kZ_ej  = 1;
constexpr double kQcorrShift = 0.208;
constexpr double kP0 = 0.618451, kP1 = 0.0080509149, kOffset = 1.022;

double omega(double x, double y, double z) {
    return std::sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z);
}
std::pair<double,double> kine_2b(double m1, double m2, double m3, double m4,
                                 double K_proj, double thetalab, double K_eject) {
    double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
    double s = m1*m1 + m2*m2 + 2*m2*Et1;
    double u = m2*m2 + m3*m3 - 2*m2*Et3;
    double m4_ex = std::sqrt((std::cos(thetalab)*omega(s,m1*m1,m2*m2)*omega(u,m2*m2,m3*m3)
                   - (s - m1*m1 - m2*m2)*(m2*m2 + m3*m3 - u))/(2*m2*m2) + s + u - m2*m2);
    double Ex = m4_ex - m4;
    double t = m2*m2 + m4_ex*m4_ex - 2*m2*Et4;
    double theta_cm = TMath::Pi() - std::acos(
        (s*s + s*(2*t - m1*m1 - m2*m2 - m3*m3 - m4_ex*m4_ex)
         + (m1*m1 - m2*m2)*(m3*m3 - m4_ex*m4_ex)) /
        (omega(s,m1*m1,m2*m2)*omega(s,m3*m3,m4_ex*m4_ex)));
    return {Ex, theta_cm*TMath::RadToDeg()};
}
}  // namespace

void fingerprint_contaminant() {
    const TString dir = "/Users/quantumlab/Downloads/C16_dp/InterpSolver_root/";
    const std::vector<TString> files = {
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
    auto H1 = [](const char* n, const char* t, int nb, double a, double b) {
        auto* h = new TH1F(n, t, nb, a, b); h->Sumw2(); return h; };
    TH1F* vz_cb=H1("vz_cb","vertex_z;vertex_z (cm);norm.",40,0,100);
    TH1F* vz_cf=H1("vz_cf","",40,0,100); TH1F* vz_rb=H1("vz_rb","",40,0,100);
    TH1F* th_cb=H1("th_cb","lab angle;#theta_{lab} (deg);norm.",45,90,180);
    TH1F* th_cf=H1("th_cf","",45,90,180); TH1F* th_rb=H1("th_rb","",45,90,180);
    TH1F* ee_cb=H1("ee_cb","proton E;E_{ej} (MeV);norm.",40,0,12);
    TH1F* ee_cf=H1("ee_cf","",40,0,12); TH1F* ee_rb=H1("ee_rb","",40,0,12);
    TH2F* ex_th=new TH2F("ex_th","E_{x} vs #theta_{cm};#theta_{cm} (deg);E_{x} (MeV)",60,0,60,80,-1,7);
    TH2F* ex_vz=new TH2F("ex_vz","E_{x} vs vertex_z;vertex_z (cm);E_{x} (MeV)",50,0,100,80,-1,7);

    long nev=0, ncb=0, nrb=0;
    for (auto& fn : files) {
        TFile* f = TFile::Open(dir+fn,"READ");
        if (!f || f->IsZombie()) { printf("skip %s\n", fn.Data()); continue; }
        TTree* T = (TTree*)f->Get("parquettree");
        if (!T) { printf("no tree %s\n", fn.Data()); f->Close(); continue; }
        double theta=0, Brho=0, zPos=0;
        T->SetBranchAddress("polar",&theta);
        T->SetBranchAddress("brho",&Brho);
        T->SetBranchAddress("vertex_z",&zPos);
        const Long64_t n = T->GetEntries();
        for (Long64_t i=0;i<n;++i) {
            T->GetEntry(i);
            const double thetaDeg = theta*TMath::RadToDeg();
            if (thetaDeg < 90.0) continue;
            const double p_ej = Brho*kZ_ej*2.99792458e2;
            const double E_ej = std::sqrt(p_ej*p_ej + kM_p*kM_p) - kM_p;
            auto pr = kine_2b(kM_C16,kM_d,kM_p,kM_C17,kEbeam,theta,E_ej);
            const double ex_raw = pr.first, theta_cm = pr.second;
            const double zCm = zPos*100.0;
            const double Ex = (ex_raw - kP1*zCm - kP0 + kOffset) + kQcorrShift;
            if (!std::isfinite(Ex) || !std::isfinite(theta_cm)) continue;
            ++nev; ex_th->Fill(theta_cm,Ex); ex_vz->Fill(zCm,Ex);
            if (theta_cm<10 || theta_cm>=40) continue;
            const bool back = theta_cm>=30.0;
            if (Ex>=0.30 && Ex<=1.10) {
                (back?vz_cb:vz_cf)->Fill(zCm); (back?th_cb:th_cf)->Fill(thetaDeg);
                (back?ee_cb:ee_cf)->Fill(E_ej); if(back) ++ncb;
            }
            if (back && Ex>=3.50 && Ex<=5.00) {
                vz_rb->Fill(zCm); th_rb->Fill(thetaDeg); ee_rb->Fill(E_ej); ++nrb;
            }
        }
        f->Close(); delete f;
    }
    printf("\n# events (polar>=90): %ld\n", nev);
    printf("# contam-window back events: %ld\n# reference back events: %ld\n", ncb, nrb);
    auto st=[](TH1F* h,const char* l){printf("  %-10s N=%6.0f mean=%7.2f rms=%6.2f\n",
             l,h->GetEntries(),h->GetMean(),h->GetRMS());};
    printf("vertex_z:\n"); st(vz_cb,"contam-bk"); st(vz_rb,"ref-bk"); st(vz_cf,"contam-fwd");
    printf("theta_lab:\n");st(th_cb,"contam-bk"); st(th_rb,"ref-bk"); st(th_cf,"contam-fwd");
    printf("E_ej:\n");     st(ee_cb,"contam-bk"); st(ee_rb,"ref-bk"); st(ee_cf,"contam-fwd");

    gStyle->SetOptStat(0);
    auto nrm=[](TH1F* h){ if(h->Integral()>0) h->Scale(1.0/h->Integral()); };
    for (TH1F* h : {vz_cb,vz_cf,vz_rb,th_cb,th_cf,th_rb,ee_cb,ee_cf,ee_rb}) nrm(h);
    TCanvas* c=new TCanvas("c","contaminant fingerprint",1500,950); c->Divide(3,2);
    auto trio=[&](int pad,TH1F* cb,TH1F* rb,TH1F* cf,const char* ti){
        c->cd(pad);
        cb->SetLineColor(kRed+1); cb->SetLineWidth(3);
        rb->SetLineColor(kAzure+1); rb->SetLineWidth(2);
        cf->SetLineColor(kGray+2); cf->SetLineWidth(2); cf->SetLineStyle(2);
        double m=std::max({cb->GetMaximum(),rb->GetMaximum(),cf->GetMaximum()});
        cb->SetMaximum(1.25*m); cb->SetTitle(ti); cb->Draw("hist");
        rb->Draw("hist same"); cf->Draw("hist same");
        auto* lg=new TLegend(0.45,0.72,0.89,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(cb,"contam, #theta_{cm}#geq30","l");
        lg->AddEntry(rb,"real 3.5-5 MeV, #geq30","l");
        lg->AddEntry(cf,"contam, #theta_{cm}<30","l"); lg->Draw();
    };
    trio(1,vz_cb,vz_rb,vz_cf,"vertex_z");
    trio(2,th_cb,th_rb,th_cf,"lab angle");
    trio(3,ee_cb,ee_rb,ee_cf,"proton energy E_{ej}");
    c->cd(4); ex_th->Draw("colz");
    c->cd(5); ex_vz->Draw("colz");
    c->cd(6); auto* tx=new TPaveText(0.04,0.05,0.96,0.95);
    tx->SetFillStyle(0); tx->SetBorderSize(0); tx->SetTextAlign(12);
    tx->AddText("0.45 MeV contaminant fingerprint");
    tx->AddText(Form("contam-back events: %ld",ncb));
    tx->AddText(Form("reference-back events: %ld",nrb));
    tx->AddText("Red vs blue: does the 0.45 MeV back-angle");
    tx->AddText("sample look like real transfer in");
    tx->AddText("vertex_z / #theta_{lab} / E_{ej}?");
    tx->AddText("Pad 4 (E_{x} vs #theta_{cm}): a slanted ridge");
    tx->AddText("crossing ~0.45 MeV at back angles = ghost.");
    tx->Draw();
    c->SaveAs("/Users/quantumlab/C17_dp_fits/plots/contam_fingerprint.png");
    printf("wrote plots/contam_fingerprint.png\n");
}
