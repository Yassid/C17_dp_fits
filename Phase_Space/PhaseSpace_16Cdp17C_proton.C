#include <TRandom.h> 

TH1F *hpstot; 
TH1F *hps1; 
TH1F *hps2; 
TH1F *hps3; 
TH1F *hps4; 
TH1F *hps5; 
TH2F *hp2d;
TRandom1 rGen(0);
TRandom1 rGen1(0);

double Ex_calculation( double ELab, double ThetaLab, double mass_in_MeV[4] , double T_b_A ) {

	double mbeam, mtarget, mejectile, mrecoil, Tbeam ;

	double MeV_to_uma = 1.0 / 931.494061 ;


	mbeam     = mass_in_MeV[0] ;
	mtarget   = mass_in_MeV[1] ;
	mejectile = mass_in_MeV[2] ;
	mrecoil   = mass_in_MeV[3] ; 

	
	Tbeam = T_b_A * mbeam * MeV_to_uma ;

	double Ebeam, pbeam, pejectile, C1, Ex ;

	Ebeam     = mbeam + Tbeam ;
	pbeam     = sqrt( 2 * mbeam * Tbeam + Tbeam * Tbeam ) ;
	pejectile = sqrt( 2 * mejectile * ELab + ELab * ELab ) ;
	C1        = mbeam * mbeam + mtarget * mtarget + mejectile * mejectile + 2 * Ebeam * mtarget ;

	Ex        = sqrt( C1 - 2 * ( Ebeam + mtarget ) * ( mejectile + ELab ) + 2 * pbeam * pejectile * cos( ThetaLab * 3.14159265359 / 180.0 )  ) - mrecoil ;

	return Ex ;
}


void Ex_and_ThetaCM_calculation_v2( double ELab, double ThetaLab, double mass_in_MeV[4], double T_b_A, double Ex_ThetaCM[2] ) {

	double ThetaCM ;

	double Eex_recoil = Ex_calculation( ELab, ThetaLab, mass_in_MeV, T_b_A ) ;

	double T_b, T3, E3, E4,  E1, Et, Et_cm, E3_cm, p1, p3_cm, p4, p3_x, p3_y, p3, theta, phi ;
	double m1, m2, m3, m4, gamma, beta, cos_theta_CM ;

	double MeV_to_uma = 1.0 / 931.494061 ;

	m1 = mass_in_MeV[0] ;
	m2 = mass_in_MeV[1] ;
	m3 = mass_in_MeV[2] ;
	m4 = mass_in_MeV[3] + Eex_recoil ; //we need ti use the one excited	

	E3 = ELab + m3 ;
	
	T_b   = T_b_A * m1 * MeV_to_uma ;
	E1    = m1 + T_b ;
	p1    = sqrt( pow( E1, 2 ) - pow( m1, 2 ) ) ;
	Et    = E1 + m2 ;
	Et_cm = sqrt( 2 * E1 * m2 + pow( m2, 2 ) + pow( m1, 2 ) ) ;
	E3_cm = 0.5 * ( Et_cm + ( pow( m3, 2 ) - pow( m4, 2 ) ) / Et_cm ) ;
	p3_cm = sqrt( pow( E3_cm, 2 ) - pow( m3, 2 ) ) ;

	gamma = Et / Et_cm ;
	beta  = p1 / Et ;

	cos_theta_CM = ( ( E3 / gamma ) - E3_cm ) / ( beta * p3_cm ) ;

	ThetaCM = 180.0 - acos( cos_theta_CM ) * 180.0 / TMath::Pi() ;

	Ex_ThetaCM[0] = Eex_recoil ; 
	Ex_ThetaCM[1] = ThetaCM    ; 
}

void PhaseSpace_1n_16Cpp16C_binWidth( int Nbin, double MinLim, double MaxLim, double sigma_Ex, double thetaCM_min, double thetaCM_max, double theta_min, double theta_max , double T_min, double T_max, int SavingDataFile, TString OutputFileName ) {

	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

	//Change masses + beam kinetic energy 
	Double_t m1,m2,m3,m4,m6,mn,bE,T1,p1,m5;

	m1 = 16.014701*931.494;	//Beam mass 10B
	m2 = 2.01410177785*931.494 ;	//Target mass d
	m3 = 1.007276 *931.494 ;	//Ejectile mass p
	m4 = 17.022586*931.494;	//Beam-like fragment 11B
	m5 = 16.014701*931.494;	// 10Be recoil
	mn = 1.008664*931.494;	//Neutron mass
	bE = 11.77;			//Beam energy MeV/A
	T1 = 16.014701*bE;		//Beam kinetic energy middle of target
	p1 = TMath::Sqrt(T1*(T1+2*m1)); //Beam momentum

	double mass_in_MeV[4] ;
	double Ex_ThetaCM[2] ;

	mass_in_MeV[0] = m1 ;
	mass_in_MeV[1] = m2 ;
	mass_in_MeV[2] = m3 ;
	mass_in_MeV[3] = m4 ;

	//px,py,pz,E in MeV
	TLorentzVector target(0.0, 0.0, 0.0, m2);
	TLorentzVector beam(0.0, 0.0, p1, T1+m1);
	TLorentzVector W = beam + target;  //LAB reference frame

	Double_t masses[3] = {m5, m3, mn};	

	//Set your Phase Space to 3 particles with masses given by array masses. 
	TGenPhaseSpace event;
	event.SetDecay(W, 3, masses);

	double etot,etot_cm,gam,beta,e3,p3,e3_cm,e3n,t3,poxy,poxyth,brho,m4ex,theta;

	etot    = T1+m1+m2;
	etot_cm = TMath::Sqrt(m1*m1+m2*m2+2*m2*(T1+m1));
	gam     = etot/etot_cm;
	beta    = TMath::Sqrt(1-1/gam/gam);

	double weight, theta_deg, ep ;

	TH1F *h_PS = new TH1F("h_PS","h_PS", Nbin, MinLim, MaxLim ) ;

	TH2F *h_kn = new TH2F("h_kn","h_kn", 180, 0, 180, 600, 0, 200 ) ;


	TFile *EfficiencySimFile ; // = new TFile("./16Cpp_phasespace/PhaseSpace_16C_pp_1n.root","recreate");
	TTree *out_tree ; 
	double ThetaCM_cal, Ex_cal, Weight_sim ; 

	if ( SavingDataFile ) {
		EfficiencySimFile = new TFile( OutputFileName,"update") ;
		out_tree = new TTree("simulated_tree", " mysim tree ") ;
		TBranch *ThetaCM_cal_b         = out_tree -> Branch( "ThetaCM_cal"        , &ThetaCM_cal        , "ThetaCM_cal/D"         ) ;
		TBranch *Ex_cal_b              = out_tree -> Branch( "Ex_cal"             , &Ex_cal             , "Ex_cal/D"              ) ;
		TBranch *Weight_sim_b          = out_tree -> Branch( "Weight_sim"         , &Weight_sim         , "Weight_sim/D"          ) ;
	}

	int I_percent = 0 ;
	int Nsim = 10000000 ;
	int factor_100 = Nsim / 100 ;
	for (Int_t n=0;  n<Nsim;n++) {


		if ( n % factor_100 == 0 ) { cout << I_percent << " /100 progress " << endl ; I_percent += 1 ; }
		//Generates a Random final state.
		weight = event.Generate();

		//Returns the Lorentz Vector for a final state.
		TLorentzVector *pOxy = event.GetDecay(0);
		TLorentzVector *pTriton= event.GetDecay(1);
		TLorentzVector *pNeut= event.GetDecay(2);

		e3n = pTriton->Energy(); //From Lorentz vector
		t3= e3n-m3 ; // MeV resolution calculated from the 3 alpha


		theta_deg = pTriton -> Theta() * 180. / TMath::Pi() ;  

		Ex_and_ThetaCM_calculation_v2( t3, theta_deg, mass_in_MeV, bE, Ex_ThetaCM ) ;

		
		
		if (theta_deg>=theta_min && theta_deg < theta_max && t3 > T_min && t3 < T_max && Ex_ThetaCM[1] < thetaCM_max && Ex_ThetaCM[1] > thetaCM_min ){ 
			double Ex_res = rGen.Gaus( Ex_ThetaCM[0] , sigma_Ex ) ;
			h_PS->Fill( Ex_res , weight );
			h_kn->Fill( theta_deg, t3 , weight );

			if ( SavingDataFile ){
				ThetaCM_cal         = Ex_ThetaCM[1] ;
				Ex_cal              = Ex_res ;
				Weight_sim          = weight ;
				out_tree->Fill() ;
			}

		}
	}

	  TCanvas *C1 = new TCanvas() ;
	  C1 -> cd() ;

	  //h_kn->Draw("colz"); 
	  h_PS->Draw("hist"); 

	  C1 -> WaitPrimitive() ;
	  C1 -> WaitPrimitive() ;

	if ( SavingDataFile ) {
		cout << "writting the data file" << endl ;
		out_tree->Write() ;
		EfficiencySimFile->Close() ;
	}


}

void PhaseSpace_2n_16Cpp16C_binWidth( int Nbin, double MinLim, double MaxLim, double sigma_Ex, double thetaCM_min, double thetaCM_max, double theta_min, double theta_max , double T_min, double T_max, int SavingDataFile, TString OutputFileName ) {

	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

	//Change masses + beam kinetic energy 
	Double_t m1,m2,m3,m4,m6,mn,bE,T1,p1,m5;

	m1 = 16.014701*931.494;	//Beam mass 15C
	m2 = 2.01410177785*931.494 ;	//Target mass p
	m3 = 1.007825*931.494 ;	//Ejectile mass p
	m4 = 17.022586*931.494;	//Beam-like fragment 15C
	m5 = 16.014701*931.494;	//Beam-like fragment 14C
	m6 = 15.010599256*931.494;	//Beam mass 14C
	mn = 1.008665*931.494;	//Neutron mass
	bE = 11.77;			//Beam energy MeV/A
	T1 = 16.014701*bE;		//Beam kinetic energy middle of target
	p1 = TMath::Sqrt(T1*(T1+2*m1)); //Beam momentum

	double mass_in_MeV[4] ;
	double Ex_ThetaCM[2] ;

	mass_in_MeV[0] = m1 ;
	mass_in_MeV[1] = m2 ;
	mass_in_MeV[2] = m3 ;
	mass_in_MeV[3] = m4 ;

	//px,py,pz,E in MeV
	TLorentzVector target(0.0, 0.0, 0.0, m2);
	TLorentzVector beam(0.0, 0.0, p1, T1+m1);
	TLorentzVector W = beam + target;  //LAB reference frame

	Double_t masses[4] = {m6, m3, mn, mn};	

	//Set your Phase Space to 3 particles with masses given by array masses. 
	TGenPhaseSpace event;
	event.SetDecay(W, 4, masses);

	double etot,etot_cm,gam,beta,e3,p3,e3_cm,e3n,t3,poxy,poxyth,brho,m4ex,theta;

	etot    = T1+m1+m2;
	etot_cm = TMath::Sqrt(m1*m1+m2*m2+2*m2*(T1+m1));
	gam     = etot/etot_cm;
	beta    = TMath::Sqrt(1-1/gam/gam);

	double weight, theta_deg, ep ;

	TH1F *h_PS = new TH1F("h_PS","h_PS", Nbin, MinLim, MaxLim ) ;

	TH2F *h_kn = new TH2F("h_kn","h_kn", 180, 0, 180, 600, 0, 200 ) ;


	TFile *EfficiencySimFile ; // = new TFile("./16Cpp_phasespace/PhaseSpace_16C_pp_1n.root","recreate");
	TTree *out_tree ; 
	double ThetaCM_cal, Ex_cal, Weight_sim ; 

	if ( SavingDataFile ) {
		EfficiencySimFile = new TFile( OutputFileName,"update") ;
		out_tree = new TTree("simulated_tree", " mysim tree ") ;
		TBranch *ThetaCM_cal_b         = out_tree -> Branch( "ThetaCM_cal"        , &ThetaCM_cal        , "ThetaCM_cal/D"         ) ;
		TBranch *Ex_cal_b              = out_tree -> Branch( "Ex_cal"             , &Ex_cal             , "Ex_cal/D"              ) ;
		TBranch *Weight_sim_b          = out_tree -> Branch( "Weight_sim"         , &Weight_sim         , "Weight_sim/D"          ) ;
	}

	int I_percent = 0 ;
	int Nsim = 10000000 ;
	int factor_100 = Nsim / 100 ;
	for (Int_t n=0;  n<Nsim;n++) {


		if ( n % factor_100 == 0 ) { cout << I_percent << " /100 progress " << endl ; I_percent += 1 ; }
		//Generates a Random final state.
		weight = event.Generate();

		//Returns the Lorentz Vector for a final state.
		TLorentzVector *pOxy    = event.GetDecay(0) ;
		TLorentzVector *pTriton = event.GetDecay(1) ;
		TLorentzVector *pNeut   = event.GetDecay(2) ;
		TLorentzVector *pNeut2  = event.GetDecay(3) ;

		e3n = pTriton->Energy(); //From Lorentz vector
		t3= e3n-m3 ; // MeV resolution calculated from the 3 alpha


		theta_deg = pTriton -> Theta() * 180. / TMath::Pi() ;  

		Ex_and_ThetaCM_calculation_v2( t3, theta_deg, mass_in_MeV, bE, Ex_ThetaCM ) ;

		
		
		if (theta_deg>=theta_min && theta_deg < theta_max && t3 > T_min && t3 < T_max && Ex_ThetaCM[1] < thetaCM_max && Ex_ThetaCM[1] > thetaCM_min ){ 
			double Ex_res = rGen.Gaus( Ex_ThetaCM[0] , sigma_Ex ) ;
			h_PS->Fill( Ex_res , weight );
			h_kn->Fill( theta_deg, t3 , weight );

			if ( SavingDataFile ){
				ThetaCM_cal         = Ex_ThetaCM[1] ;
				Ex_cal              = Ex_res ;
				Weight_sim          = weight ;
				out_tree->Fill() ;
			}

		}
	}

	  TCanvas *C1 = new TCanvas() ;
	  C1 -> cd() ;

	  //h_kn->Draw("colz"); 
	  h_PS->Draw("hist"); 

	  std::cout << "Finished Drawing Phase Space 2n" << std::endl ;

	  C1 -> WaitPrimitive() ;
	  C1 -> WaitPrimitive() ;

	  std::cout << "Finished Drawing Phase Space 2n" << std::endl ;

	if ( SavingDataFile ) {
		cout << "writting the data file" << endl ;
		out_tree->Write() ;
		EfficiencySimFile->Close() ;
		std::cout << "Finished Writing Phase Space 2n Data File" << std::endl ;
	}


}

void PhaseSpace_16Cdp17C_proton() {
	double MinLim, MaxLim ;
	int Nbin ;

	MinLim =  -1. ;
	MaxLim = 24. ; 
	Nbin   = 250 ;

	//Acceptance and cuts 
	double theta_min  = 90.0 ;
	double theta_max  = 145.0 ;

	double thetaCM_min  = 10.0 ;
	double thetaCM_max  = 65.0 ;

	double T_min  = 0.5 ;
	double T_max  = 100.0 ;

	double sigma_Ex = 0.200 ; //THis can be modified in the function if the resolution dependes on the Ex energy (or directly include the E anf angle resolution if preferred)

	int SavingDataFile = 1 ;
	TString OutputFileName_1n = "./PhaseSpace_16C_dp_1n.root" ;
	PhaseSpace_1n_16Cpp16C_binWidth( Nbin, MinLim, MaxLim, sigma_Ex, thetaCM_min, thetaCM_max, theta_min, theta_max, T_min, T_max, SavingDataFile, OutputFileName_1n  ) ;
	TString OutputFileName_2n = "./PhaseSpace_16C_dp_2n.root" ;
	PhaseSpace_2n_16Cpp16C_binWidth( Nbin, MinLim, MaxLim, sigma_Ex, thetaCM_min, thetaCM_max, theta_min, theta_max, T_min, T_max, SavingDataFile, OutputFileName_2n  ) ;

}
