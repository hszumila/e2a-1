/*
// ===============================================================================
Code to calculate parameters needed for neutron #beta corrections

We take as input 3He(e,e’pp) events
(i.e. # particles in final state >=3 and # protons in final state >= 2)

we only put these events through the following cuts before using them here:
|	Ignore events that have no particle candidates
|	Electron ID:
|	DC, CC, EC, SC status good, q < 0
|	Cut on E/p vs p, cut on min and max momentum, Energy deposited in EC inner layer, Total energy deposited in EC
|	Momentum of electron candidate
|	Only for 2 GeV data:
|		Momentum of electron candidate < Eb
|		Delta t SC - CC (sector dependent) > some value
|		angle between CC hit and nearest SC cut < 0.1 rad
|	Electron corrections: momentum corrections, vz corrections
|	
|	Positive ID: q > 0, SC status, DC status, Global status good
|	
|	Proton ID:id = 2212, Cut on Delta T vs p distribution
|	Proton Corrections: energy loss correction, vz corrections
|
|	Neutral ID: global status is good. Beta < 1
|
|	No fiducial cuts applied
 */

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include <iostream>

using namespace std;

const double Me      = .000511;         // Electron     mass in GeV
const double Mp      = 0.93827;         // Proton       mass in GeV
const double Mn      = 0.93957;         // Neutron      mass in GeV
const double Mpip    = 0.13957;         // Charged pion mass in GeV
const double He3_be  = 0.00770;         // 3He binding energy in GeV
const double M_tar   = 2*Mp+Mn-He3_be;  // 3Helium  mass in GeV (target)
const double c_cm_s  = 2.99792E+10; //speed of the light in cm/s
const double ns2s    = 1.0E-09 ;

double max_outta_three ( double A , double B, double C );
Double_t func_fit(Double_t *x, Double_t *par);
void fit_gfunction(TH1F * gPlot, float minx, float maxx);

// ====================================================================================================================================================
// Main Function
int main() {

	//  option = 0: 4.4 GeV data
	//         = 1: 2.2 GeV data
	int option = 0;

	bool eppn, eppnpipi;

	//Define Branch variables and other variables needed to determine parameters
	int nEntries;

	const double Mn_up   = 1.05000;		// Delta neutron mass in GeV for mass distribution cut
	const double Mn_do   = 0.85000;		// Delta neutron mass in GeV for mass distribution cut

	double Ebeam, M_miss_eppn, M_miss_eppnpipi, nBeta, Rcalc, cos_th_miss_ec_eppn, cos_th_miss_ec_eppnpipi;
	TVector3 EC_hit_pos, P_miss_eppn, P_miss_eppnpipi, u1;
	TLorentzVector mom_e,mom_p1,mom_p2, E0, V4_tar, mom_pip, mom_pim;
	double sig_e_p1, sig_e_p2, sig_p1p2, sig_e_pip, sig_e_pim, sig_pip_pim, temp1, temp2,
	       Z_e_p1_abs, Z_e_p2_abs, Z_p1p2_abs, Z_e_pip_abs, Z_e_pim_abs, Z_pip_pim_abs, vz_o_sig;
	// --------------------------------------------------------------------------
	//Read in data files	
	TChain *t = new TChain("T");
	if(option==0){
		// 4.4 GeV parameters
		cout << "Creating file for 4.4 GeV data" << endl;

		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18443.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18449.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18453.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18454.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18458.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18459.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18460.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18461.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18462.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18463.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18464.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18465.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18469.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18470.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18471.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18472.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18473.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18474.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18475.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18476.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18477.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18520.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18521.root");

		Ebeam = 4.4; //Beam energy in GeV
	}
	else if(option==1){
		// 2.2 GeV parameters
		cout << "Creating file for 2.2 GeV data" << endl;

		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18346.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18349.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18350.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18351.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18365.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18366.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18367.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18368.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18370.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18374.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18380.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18381.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18382.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18383.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18384.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18385.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18386.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18388.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18389.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18402.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18405.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18406.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18407.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18408.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18409.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18410.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18415.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18416.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18428.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18430.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18431.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18432.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18433.root");
		t->Add("../3_Neutron_path_length_corr/data/skimmed_data_18437.root");

		Ebeam = 2.2; //Beam energy in GeV
	}
	else{
		cout << "Wrong option. Check variable 'option'." << endl;
		exit(1);
	}
	// --------------------------------------------------------------------------
	int nParticles, nProtons, nNeutrons;
	int Part_type [20], stat_dc  [20], stat_sc[20], stat_ec[20];
	double t0;
	double mom_x  [20], mom_y    [20], mom_z  [20], Mass   [20], charge[20], beta[20], vtx_z_cor[20],
	       ec_time[20], ec_path  [20], ec_in  [20], ec_out [20], ec_tot[20], ec_x[20], ec_y     [20], ec_z[20],
	       ec_u   [20], ec_v     [20], ec_w   [20];
	//Define Branch Addresses
	t->SetBranchAddress("nParticles",&nParticles);
	t->SetBranchAddress("nProtons"  ,&nProtons  );
	t->SetBranchAddress("nNeutrons" ,&nNeutrons );
	t->SetBranchAddress("t0"        ,&t0        );
	t->SetBranchAddress("charge"    , charge    );
	t->SetBranchAddress("Part_type" , Part_type );
	t->SetBranchAddress("stat_dc"   , stat_dc   );
	t->SetBranchAddress("stat_sc"   , stat_sc   );
	t->SetBranchAddress("stat_ec"   , stat_ec   );
	t->SetBranchAddress("ec_time"	, ec_time   );
	t->SetBranchAddress("ec_path"	, ec_path   );
	t->SetBranchAddress("ec_x"      , ec_x      );
	t->SetBranchAddress("ec_y"      , ec_y      );
	t->SetBranchAddress("ec_z"      , ec_z      );
	t->SetBranchAddress("ec_u"      , ec_u      );
	t->SetBranchAddress("ec_v"      , ec_v      );
	t->SetBranchAddress("ec_w"      , ec_w      );
	t->SetBranchAddress("ec_in"	, ec_in     );
	t->SetBranchAddress("ec_out"	, ec_out    );
	t->SetBranchAddress("ec_tot"	, ec_tot    );
	t->SetBranchAddress("mom_x"     , mom_x     );
	t->SetBranchAddress("mom_y"     , mom_y     );
	t->SetBranchAddress("mom_z"     , mom_z     );
	t->SetBranchAddress("Mass"      , Mass      );
	t->SetBranchAddress("beta"      , beta      );
	t->SetBranchAddress("vtx_z_cor" , vtx_z_cor );

	// --------------------------------------------------------------------------
	// Define Histograms
	// -----------------------------------
	// EC coordinates cut
	TH1D * h1_u_eppn        = new TH1D("h1_u_eppn"     ,"^{3}He(e,e'ppn);EC_{u} [cm];Counts"                 ,100,   0., 500.);
	TH1D * h1_v_eppn        = new TH1D("h1_v_eppn"     ,"^{3}He(e,e'ppn);EC_{v} [cm];Counts"                 ,100,   0., 500.);
	TH1D * h1_w_eppn        = new TH1D("h1_w_eppn"     ,"^{3}He(e,e'ppn);EC_{w} [cm];Counts"                 ,100,   0., 500.);
	TH1D * h1_u_eppn_cut    = new TH1D("h1_u_eppn_cut" ,"^{3}He(e,e'ppn);EC_{u} [cm];Counts"                 ,100,   0., 500.);
	TH1D * h1_v_eppn_cut    = new TH1D("h1_v_eppn_cut" ,"^{3}He(e,e'ppn);EC_{v} [cm];Counts"                 ,100,   0., 500.);
	TH1D * h1_w_eppn_cut    = new TH1D("h1_w_eppn_cut" ,"^{3}He(e,e'ppn);EC_{w} [cm];Counts"                 ,100,   0., 500.);
	h1_u_eppn_cut -> SetFillColor(2);
	h1_u_eppn_cut -> SetLineColor(2);
	h1_v_eppn_cut -> SetFillColor(2);
	h1_v_eppn_cut -> SetLineColor(2);
	h1_w_eppn_cut -> SetFillColor(2);
	h1_w_eppn_cut -> SetLineColor(2);
	// ---
	TH1D * h1_u_eppnpipi     = new TH1D("h1_u_eppnpipi"     ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{u} [cm];Counts"  ,100,   0., 500.);
	TH1D * h1_v_eppnpipi     = new TH1D("h1_v_eppnpipi"     ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{v} [cm];Counts"  ,100,   0., 500.);
	TH1D * h1_w_eppnpipi     = new TH1D("h1_w_eppnpipi"     ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{w} [cm];Counts"  ,100,   0., 500.);
	TH1D * h1_u_eppnpipi_cut = new TH1D("h1_u_eppnpipi_cut" ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{u} [cm];Counts"  ,100,   0., 500.);
	TH1D * h1_v_eppnpipi_cut = new TH1D("h1_v_eppnpipi_cut" ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{v} [cm];Counts"  ,100,   0., 500.);
	TH1D * h1_w_eppnpipi_cut = new TH1D("h1_w_eppnpipi_cut" ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});EC_{w} [cm];Counts"  ,100,   0., 500.);
	h1_u_eppnpipi_cut -> SetFillColor(2);
	h1_u_eppnpipi_cut -> SetLineColor(2);
	h1_v_eppnpipi_cut -> SetFillColor(2);
	h1_v_eppnpipi_cut -> SetLineColor(2);
	h1_w_eppnpipi_cut -> SetFillColor(2);
	h1_w_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------
	// Missing mass cut
	TH1F *h1_M_miss_eppn         = new TH1F("h1_M_miss_eppn"         ,";m_{miss} [GeV];Counts",100,   0., 2.5);
	TH1F *h1_M_miss_eppn_cut     = new TH1F("h1_M_miss_eppn_cut"     ,""                      ,100,   0., 2.5);
	h1_M_miss_eppn_cut -> SetFillColor(2);
	h1_M_miss_eppn_cut -> SetLineColor(2);
	// ---
	TH1F *h1_M_miss_eppnpipi     = new TH1F("h1_M_miss_eppnpipi"     ,";m_{miss} [GeV];Counts",100,   0., 2.5);
	TH1F *h1_M_miss_eppnpipi_cut = new TH1F("h1_M_miss_eppnpipi_cut" ,""                      ,100,   0., 2.5);
	h1_M_miss_eppnpipi_cut -> SetFillColor(2);
	h1_M_miss_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------
	// Pmiss cut
	TH1F *h1_P_miss_eppn         = new TH1F("h1_P_miss_eppn"         ,";p_{miss} [GeV];Counts", 50,   0., 2.5);
	TH1F *h1_P_miss_eppn_cut     = new TH1F("h1_P_miss_eppn_cut"     ,""                      , 50,   0., 2.5);
	h1_P_miss_eppn_cut -> SetFillColor(2);
	h1_P_miss_eppn_cut -> SetLineColor(2);
	// ---
	TH1F *h1_P_miss_eppnpipi     = new TH1F("h1_P_miss_eppnpipi"     ,";p_{miss} [GeV];Counts", 50,   0., 2.5);
	TH1F *h1_P_miss_eppnpipi_cut = new TH1F("h1_P_miss_eppnpipi_cut" ,""                      , 50,   0., 2.5);
	h1_P_miss_eppnpipi_cut -> SetFillColor(2);
	h1_P_miss_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------
	// Vertex difference cut
	TH1F *h1_delta_vz_e_p1    = new TH1F("h1_delta_vz_e_p1"   ,"#Delta v_{Z};e - p1;Counts"     , 50,  -5., 5. );
	TH1F *h1_delta_vz_e_p2    = new TH1F("h1_delta_vz_e_p2"   ,"#Delta v_{Z};e - p2;Counts"     , 50,  -5., 5. );
	TH1F *h1_delta_vz_p1_p2   = new TH1F("h1_delta_vz_p1_p2"  ,"#Delta v_{Z};p1 - p2;Counts"    , 50,  -5., 5. );
	TH1F *h1_delta_vz_e_pip   = new TH1F("h1_delta_vz_e_pip"  ,"#Delta v_{Z};e - pi+;Counts"    , 50,  -5., 5. );
	TH1F *h1_delta_vz_e_pim   = new TH1F("h1_delta_vz_e_pim"  ,"#Delta v_{Z};e - pi-;Counts"    , 50,  -5., 5. );
	TH1F *h1_delta_vz_pip_pim = new TH1F("h1_delta_vz_pip_pim","#Delta v_{Z};pi+ - pi-;Counts"  , 50,  -5., 5. );
	// ---
	TH1F *h1_max_vz_sig_eppn     = new TH1F("h1_max_vz_sig_eppn"    ,";max(|Z_{e}-Z_{p1}|,|Z_{e}-Z_{p2}|,|Z_{p1}-Z_{p2}|)/#sigma;Counts",50,0.,6.);
	TH1F *h1_max_vz_sig_eppn_cut = new TH1F("h1_max_vz_sig_eppn_cut",";max(|Z_{e}-Z_{p1}|,|Z_{e}-Z_{p2}|,|Z_{p1}-Z_{p2}|)/#sigma;Counts",50,0.,6.);
	h1_max_vz_sig_eppn_cut -> SetFillColor(2);
	h1_max_vz_sig_eppn_cut -> SetLineColor(2);
	// ---
	TH1F *h1_max_vz_sig_eppnpipi     = new TH1F("h1_max_vz_sig_eppnpipi"    ,";max(|Z_{e}-Z_{p1}|,|Z_{e}-Z_{p2}|,|Z_{p1}-Z_{p2}|)/#sigma;Counts",50,0.,6.);
	TH1F *h1_max_vz_sig_eppnpipi_cut = new TH1F("h1_max_vz_sig_eppnpipi_cut",";max(|Z_{e}-Z_{p1}|,|Z_{e}-Z_{p2}|,|Z_{p1}-Z_{p2}|)/#sigma;Counts",50,0.,6.);
	h1_max_vz_sig_eppnpipi_cut -> SetFillColor(2);
	h1_max_vz_sig_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------
	// Beta cut
	TH1F *h1_beta_eppn      = new TH1F("h1_beta_eppn"      ,";#beta;Counts", 50,   0., 1.0);
	TH1F *h1_beta_eppn_cut  = new TH1F("h1_beta_eppn_cut"  ,";#beta;Counts", 50,   0., 1.0);
	h1_beta_eppn_cut -> SetFillColor(2);
	h1_beta_eppn_cut -> SetLineColor(2);
	// ---
	TH1F *h1_beta_eppnpipi     = new TH1F("h1_beta_eppnpipi"     ,";#beta;Counts", 50,   0., 1.0);
	TH1F *h1_beta_eppnpipi_cut = new TH1F("h1_beta_eppnpipi_cut" ,";#beta;Counts", 50,   0., 1.0);
	h1_beta_eppnpipi_cut -> SetFillColor(2);
	h1_beta_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------
	// cos( miss , ec ) cut
	TH1F *h1_cos_miss_ec_eppn    = new TH1F("h1_cos_miss_ec_eppn"    ,"^{3}He(e,e'ppn);Cos#theta_{miss,EC};Counts",100, 0.95,1.001);
	TH1F *h1_cos_miss_ec_eppn_cut= new TH1F("h1_cos_miss_ec_eppn_cut",""                                          ,100, 0.95,1.001);
	h1_cos_miss_ec_eppn_cut -> SetFillColor(2);
	h1_cos_miss_ec_eppn_cut -> SetLineColor(2);
	// ---
	TH1F *h1_cos_miss_ec_eppnpipi     = new TH1F("h1_cos_miss_ec_eppnpipi"    ,"^{3}He(e,e'ppn#pi^{+}#pi^{-});Cos#theta_{miss,EC};Counts",100, 0.95,1.001);
	TH1F *h1_cos_miss_ec_eppnpipi_cut = new TH1F("h1_cos_miss_ec_eppnpipi_cut",""                                                        ,100, 0.95,1.001);
	h1_cos_miss_ec_eppnpipi_cut -> SetFillColor(2);
	h1_cos_miss_ec_eppnpipi_cut -> SetLineColor(2);
	// -----------------------------------

	TH1F * h1_pmiss_epp        = new TH1F("h1_pmiss_epp"        , "^{3}He(e,e'pp)n;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
        TH1F * h1_pmiss_epppipi    = new TH1F("h1_pmiss_epppipi"    , "^{3}He(e,e'pp#pi^{+}#pi^{-})n;p_{miss} [GeV];Counts"   , 20 , 0.5 , 2.1 );
	TH1F * h1_pmiss_both_inc   = new TH1F("h1_pmiss_both_inc"   , "both, inclusive;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );

	TH1F * h1_pmiss_eppn       = new TH1F("h1_pmiss_eppn"       , "^{3}He(e,e'ppn);p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
	TH1F * h1_pmiss_epppipin   = new TH1F("h1_pmiss_epppipin"   , "^{3}He(e,e'pp#pi^{+}#pi^{-}n);p_{miss} [GeV];Counts"   , 20 , 0.5 , 2.1 );
	TH1F * h1_pmiss_both_exc   = new TH1F("h1_pmiss_both_exc"   , "both, exclusive;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );

	TH1F * h1_det_eff_eppn     = new TH1F("h1_det_eff_eppn"     , ""                                                      , 20 , 0.5 , 2.1 );
	TH1F * h1_det_eff_eppnpipi = new TH1F("h1_det_eff_eppnpipi" , ""                                                      , 20 , 0.5 , 2.1 );
	TH1F * h1_det_eff_both     = new TH1F("h1_det_eff_both"     , ""                                                      , 20 , 0.5 , 2.1 );

	cout << "Histograms have been defined" << endl;

	// ==============================================================================================================================================
	E0.SetXYZT    (0.,0.,sqrt(Ebeam*Ebeam-Me*Me),Ebeam);
	V4_tar.SetXYZT(0.,0.,0.                     ,M_tar);

	nEntries = t->GetEntries();

	cout << "Looping over events the first time" << endl;
	for(int evt=0;evt<nEntries;evt++){
		if(evt%100000==0) cout << " events processed = " << evt << ", out of " << nEntries << endl;
		t->GetEntry(evt);

		// --------------------------------------------------------------------------
		// Make sure the input root file has the correct particle content       
		if(nParticles<3){
			cout << "Error: particles in the input file; this code should use 3He(e,e'pp) as input." << endl;
			exit(1);
		}

		int p1 = 0;
		int p2 = 0;
		int pip= 0;
		int pim= 0;
		int n1 = 0;

		for(int i = 1 ; i < nParticles ; i++){
			if (Part_type[i] == 2212){
				if       (p1==0)                p1 = i;
				else if ((p1!=0)&&(p2==0))      p2 = i;
			}
			// Check if there's a pion+ in the event
			if ((Part_type[i] ==  211)&&(pip==0)){
				pip = i;
			}
			// Check if there's a pion- in the event
			if ((Part_type[i] == -211)&&(pim==0)){
				pim = i;
			}

			// neutral particles with good conditions, not just for neutron.
			else if ((charge[i] == 0)&&(n1==0)){
				n1 = i; 
			}
		}

		if (p1==0||p2==0){cout << "No two protons in this event" << endl;       continue;}

		h1_delta_vz_e_p1  -> Fill(vtx_z_cor[ 0]-vtx_z_cor[p1]);
		h1_delta_vz_e_p2  -> Fill(vtx_z_cor[ 0]-vtx_z_cor[p2]);
		h1_delta_vz_p1_p2 -> Fill(vtx_z_cor[p1]-vtx_z_cor[p2]);

		if((pip!=0)&&(pim!=0)){
			h1_delta_vz_e_pip   -> Fill(vtx_z_cor[  0]-vtx_z_cor[pip]);
			h1_delta_vz_e_pim   -> Fill(vtx_z_cor[  0]-vtx_z_cor[pim]);
			h1_delta_vz_pip_pim -> Fill(vtx_z_cor[pip]-vtx_z_cor[pim]);
		}

	}
	// --------------------------------------------------
	TCanvas * c0 = new TCanvas("c0","#Delta vz",1200,900);
	c0 -> Divide(3,2);
	c0 -> cd(1);    h1_delta_vz_e_p1    -> Draw();
	c0 -> cd(2);    h1_delta_vz_e_p2    -> Draw();
	c0 -> cd(3);    h1_delta_vz_p1_p2   -> Draw();
	c0 -> cd(4);	h1_delta_vz_e_pip   -> Draw();
	c0 -> cd(5);	h1_delta_vz_e_pim   -> Draw();
	c0 -> cd(6);	h1_delta_vz_pip_pim -> Draw();

	TFitResultPtr f_vz_e_p1    = h1_delta_vz_e_p1    -> Fit("gaus","S", "", -3., 3.);
	TFitResultPtr f_vz_e_p2    = h1_delta_vz_e_p2    -> Fit("gaus","S", "", -3., 3.);
	TFitResultPtr f_vz_p1p2    = h1_delta_vz_p1_p2   -> Fit("gaus","S", "", -3., 3.);
	TFitResultPtr f_vz_e_pip   = h1_delta_vz_e_pip   -> Fit("gaus","S", "", -3., 3.);
	TFitResultPtr f_vz_e_pim   = h1_delta_vz_e_pim   -> Fit("gaus","S", "", -3., 3.);
	TFitResultPtr f_vz_pip_pim = h1_delta_vz_pip_pim -> Fit("gaus","S", "", -3., 3.);

	sig_e_p1    = f_vz_e_p1    -> Value(2);
	sig_e_p2    = f_vz_e_p2    -> Value(2);
	sig_p1p2    = f_vz_p1p2    -> Value(2);
	sig_e_pip   = f_vz_e_pip   -> Value(2);
	sig_e_pim   = f_vz_e_pim   -> Value(2);
	sig_pip_pim = f_vz_pip_pim -> Value(2);

	// ======================================================================================================
	cout << "Looping over events the second time" << endl;
	for(int evt=0;evt<nEntries;evt++){

		if(evt%100000==0) cout << " events processed = " << evt << ", out of " << nEntries << endl;
		t->GetEntry(evt);

		int p1 = 0;
		int p2 = 0;
		int n1 = 0;

		for(int i = 1 ; i < nParticles ; i++){
			if (Part_type[i] == 2212){
				if       (p1==0)		p1 = i;
				else if ((p1!=0)&&(p2==0))	p2 = i;
			}
			else if ((charge[i] == 0)&&(n1==0))	n1 = i;
		}

		if (p1==0||p2==0){cout << "No two protons in this event" << endl;	continue;}

		// --------------------------------------------------------------------------
		mom_e .SetXYZM(mom_x[ 0],mom_y[ 0],mom_z[ 0],Me);               // Electron
		mom_p1.SetXYZM(mom_x[p1],mom_y[p1],mom_z[p1],Mp);               // Proton 1
		mom_p2.SetXYZM(mom_x[p2],mom_y[p2],mom_z[p2],Mp);               // Proton 2

		// Calculate missing mass and momentum	
		M_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2).M   ();
		P_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2).Vect();

		// --------------------------------------------------------------------------
		// Cuts
		// -----------------------------------
		// 3He(e,e'ppn)
		// -----------------------------------
		// Missing mass cut	
		h1_M_miss_eppn     -> Fill(M_miss_eppn);
		if((M_miss_eppn < Mn_up)&&(M_miss_eppn > Mn_do)){
			h1_M_miss_eppn_cut -> Fill(M_miss_eppn);

			// -----------------------------------
			// Vertex difference cut
			Z_e_p1_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p1])/sig_e_p1;
			Z_e_p2_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p2])/sig_e_p2;
			Z_p1p2_abs = TMath::Abs(vtx_z_cor[p1]-vtx_z_cor[p2])/sig_p1p2;
			vz_o_sig = max_outta_three(Z_e_p1_abs,Z_e_p2_abs,Z_p1p2_abs);
			h1_max_vz_sig_eppn     -> Fill(vz_o_sig);
			if(vz_o_sig<3.){
				h1_max_vz_sig_eppn_cut -> Fill(vz_o_sig);	

				// -----------------------------------
				// Pmiss cut
				h1_P_miss_eppn     -> Fill(P_miss_eppn.Mag());
				if(P_miss_eppn.Mag() > 0.5) {
					h1_P_miss_eppn_cut -> Fill(P_miss_eppn.Mag());

					h1_pmiss_epp        -> Fill(P_miss_eppn    .Mag());
					h1_pmiss_both_inc   -> Fill(P_miss_eppn    .Mag());

					// EXCLUSIVE REACTION BEYOND THIS POINT

					if (n1!=0){ //No neutral hit in the EC

						if( (stat_ec[n1]>0) && (stat_dc[n1]<=0) && (stat_sc[n1]<=0) ){

							// Cut 10 cm from edges of calorimeter
							h1_u_eppn -> Fill(ec_u[n1]);
							h1_v_eppn -> Fill(ec_v[n1]);
							h1_w_eppn -> Fill(ec_w[n1]);
							if ( (ec_u[n1] > 40) && (ec_v[n1] < 360 ) && (ec_w[n1] < 395) ){
								h1_u_eppn_cut -> Fill(ec_u[n1]);
								h1_v_eppn_cut -> Fill(ec_v[n1]);
								h1_w_eppn_cut -> Fill(ec_w[n1]);
								// -----------------------------------
								// Beta cut
								h1_beta_eppn     -> Fill(beta[n1]);
								if(beta[n1] < 0.95){
									h1_beta_eppn_cut -> Fill(beta[n1]);
									// -----------------------------------
									// cos( miss , ec ) cut
									EC_hit_pos.SetXYZ(ec_x[n1],ec_y[n1],ec_z[n1]);
									cos_th_miss_ec_eppn = P_miss_eppn*EC_hit_pos/P_miss_eppn.Mag()/EC_hit_pos.Mag();	
									h1_cos_miss_ec_eppn     -> Fill(cos_th_miss_ec_eppn);
									if(cos_th_miss_ec_eppn > 0.995) {
										h1_cos_miss_ec_eppn_cut -> Fill(cos_th_miss_ec_eppn);
										// --------------------------------------------------------------------------
										// Fill histograms
										if((ec_tot[n1]>0)&&(ec_time[n1]>0)){
											h1_pmiss_eppn     -> Fill(P_miss_eppn.Mag());
        										h1_pmiss_both_exc -> Fill(P_miss_eppn.Mag());
											h1_det_eff_eppn   -> Fill(P_miss_eppn.Mag());
                                                                                        h1_det_eff_both   -> Fill(P_miss_eppn.Mag());
										}
									}
								}
							}
						}
					}
				}
			}
		}	
	}
	cout << "Done with the loop" << endl;

	// -----------------------------------
	// 3He(e,e'ppn pi+pi-)
	// -----------------------------------
	// ======================================================================================================
	cout << "Looping over events the third time" << endl;
	for(int evt=0;evt<nEntries;evt++){

		eppnpipi = false;

		if(evt%100000==0) cout << " events processed = " << evt << ", out of " << nEntries << endl;
		t->GetEntry(evt);

		int p1 = 0;
		int p2 = 0;
		int pip= 0;
		int pim= 0;
		int n1 = 0;

		for(int i = 1 ; i < nParticles ; i++){
			if (Part_type[i] == 2212){
				if       (p1==0)                p1 = i;
				else if ((p1!=0)&&(p2==0))      p2 = i;
			}
			// Check if there's a pion+ in the event
			if ((Part_type[i] ==  211)&&(pip==0))	pip = i;
			// Check if there's a pion- in the event
			if ((Part_type[i] == -211)&&(pim==0))	pim = i;

			// neutral particles with good conditions, not just for neutron.
			else if ((charge[i] == 0)&&(n1==0))	n1 = i;
		}

		if (  p1==0 ||  p2 ==0 ){cout << "No two protons in this event" << endl;       continue;}
		if ( pip==0 || pim ==0 ) continue;

		if((p1!=0)&&(p2!=0)&&(pip!=0)&&(pim!=0)){eppnpipi = true;}
		// --------------------------------------------------------------------------
		mom_e .SetXYZM(mom_x[ 0],mom_y[ 0],mom_z[ 0],Me);               // Electron
		mom_p1.SetXYZM(mom_x[p1],mom_y[p1],mom_z[p1],Mp);               // Proton 1
		mom_p2.SetXYZM(mom_x[p2],mom_y[p2],mom_z[p2],Mp);               // Proton 2

		// Calculate missing mass and momentum
		if(eppnpipi){
			mom_pip.SetXYZM(mom_x[pip],mom_y[pip],mom_z[pip],Mpip);               // Pi+
			mom_pim.SetXYZM(mom_x[pim],mom_y[pim],mom_z[pim],Mpip);               // Pi-

			M_miss_eppnpipi = (E0 + V4_tar - mom_e - mom_p1 - mom_p2 - mom_pip - mom_pim).M   ();
			P_miss_eppnpipi = (E0 + V4_tar - mom_e - mom_p1 - mom_p2 - mom_pip - mom_pim).Vect();
		}

		// --------------------------------------------------------------------------
		// Cuts
		// -----------------------------------
		// Missing mass cut     
		h1_M_miss_eppnpipi     -> Fill(M_miss_eppnpipi);
		if((M_miss_eppnpipi < Mn_up)&&(M_miss_eppnpipi > Mn_do)){
			h1_M_miss_eppnpipi_cut -> Fill(M_miss_eppnpipi);

			// -----------------------------------
			// Vertex difference cut
			Z_e_p1_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p1])/sig_e_p1   ;
			Z_e_p2_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p2])/sig_e_p2   ;
			Z_p1p2_abs    = TMath::Abs(vtx_z_cor[ p1]-vtx_z_cor[ p2])/sig_p1p2   ;
			Z_e_pip_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pip])/sig_e_pip  ;
			Z_e_pim_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pim])/sig_e_pim  ;
			Z_pip_pim_abs = TMath::Abs(vtx_z_cor[pip]-vtx_z_cor[pim])/sig_pip_pim;
			temp1    = max_outta_three(Z_e_p1_abs ,Z_e_p2_abs ,Z_p1p2_abs   );
			temp2    = max_outta_three(Z_e_pip_abs,Z_e_pim_abs,Z_pip_pim_abs);
			vz_o_sig = TMath::Max(temp1,temp2);
			h1_max_vz_sig_eppnpipi     -> Fill(vz_o_sig);
			if(vz_o_sig<3.){
				h1_max_vz_sig_eppnpipi_cut -> Fill(vz_o_sig);

				// -----------------------------------
				// Pmiss cut
				h1_P_miss_eppnpipi     -> Fill(P_miss_eppnpipi.Mag());
				if(P_miss_eppnpipi.Mag() > 0.5){
					h1_P_miss_eppnpipi_cut -> Fill(P_miss_eppnpipi.Mag());

					h1_pmiss_epppipi    -> Fill(P_miss_eppnpipi.Mag());
        				h1_pmiss_both_inc   -> Fill(P_miss_eppnpipi.Mag());

					// EXCLUSIVE REACTION BEYOND THIS POINT

					if (n1!=0){ //No neutral hit in the EC

						if( (stat_ec[n1]>0) && (stat_dc[n1]<=0) && (stat_sc[n1]<=0) ){

							// Cut 10 cm from edges of calorimeter
							h1_u_eppnpipi -> Fill(ec_u[n1]);
							h1_v_eppnpipi -> Fill(ec_v[n1]);
							h1_w_eppnpipi -> Fill(ec_w[n1]);
							if ( (ec_u[n1] > 40) || (ec_v[n1] < 360 ) || (ec_w[n1] < 395) ){
								h1_u_eppnpipi_cut -> Fill(ec_u[n1]);
								h1_v_eppnpipi_cut -> Fill(ec_v[n1]);
								h1_w_eppnpipi_cut -> Fill(ec_w[n1]);
								// -----------------------------------
								// Beta cut
								h1_beta_eppnpipi     -> Fill(beta[n1]);
								if(beta[n1] < 0.95){
									h1_beta_eppnpipi_cut -> Fill(beta[n1]);
									// -----------------------------------
									// cos( miss , ec ) cut
									EC_hit_pos.SetXYZ(ec_x[n1],ec_y[n1],ec_z[n1]);
									cos_th_miss_ec_eppnpipi = P_miss_eppnpipi*EC_hit_pos/P_miss_eppnpipi.Mag()/EC_hit_pos.Mag();
									h1_cos_miss_ec_eppnpipi     -> Fill(cos_th_miss_ec_eppnpipi);
									if(cos_th_miss_ec_eppnpipi > 0.995) {
										h1_cos_miss_ec_eppnpipi_cut -> Fill(cos_th_miss_ec_eppnpipi);
										// --------------------------------------------------------------------------
										// Fill histograms
										if((ec_tot[n1]>0)&&(ec_time[n1]>0)){
											h1_pmiss_epppipin   -> Fill(P_miss_eppnpipi.Mag());
        										h1_pmiss_both_exc   -> Fill(P_miss_eppnpipi.Mag());
											h1_det_eff_eppnpipi -> Fill(P_miss_eppnpipi.Mag());
											h1_det_eff_both     -> Fill(P_miss_eppnpipi.Mag());
										}
									}
								}
							}	
						}
					}
				}
			}
		}
		// --------------------------------------------------------------------------


	}
	cout << "Done with the loop" << endl;

	// --------------------------------------------------------------------------
	// Calculating neutron detection efficiency
	h1_det_eff_eppn     -> Divide( h1_pmiss_epp      );
        h1_det_eff_eppnpipi -> Divide( h1_pmiss_epppipi  );
        h1_det_eff_both     -> Divide( h1_pmiss_both_inc );

	double error_eppn, error_eppnpipi, error_both, A, B, eA, eB;
	
	for( int i = 1 ; i <= h1_det_eff_eppn -> GetXaxis()->GetNbins(); i++){
		A  = h1_pmiss_eppn -> GetBinContent(i);
		B  = h1_pmiss_epp  -> GetBinContent(i);
		eA = h1_pmiss_eppn -> GetBinError  (i);
		eB = h1_pmiss_epp  -> GetBinError  (i);
		error_eppn = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
		h1_det_eff_eppn -> SetBinError(i,error_eppn);
	}

	for( int i = 1 ; i <= h1_det_eff_eppnpipi -> GetXaxis()->GetNbins(); i++){
                A  = h1_pmiss_epppipin -> GetBinContent(i);
                B  = h1_pmiss_epppipi  -> GetBinContent(i);
                eA = h1_pmiss_epppipin -> GetBinError  (i);
                eB = h1_pmiss_epppipi  -> GetBinError  (i);
                error_eppnpipi = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
                h1_det_eff_eppnpipi -> SetBinError(i,error_eppnpipi);
        }

	for( int i = 1 ; i <= h1_det_eff_both -> GetXaxis()->GetNbins(); i++){
                A  = h1_pmiss_both_exc -> GetBinContent(i);
                B  = h1_pmiss_both_inc -> GetBinContent(i);
                eA = h1_pmiss_both_exc -> GetBinError  (i);
                eB = h1_pmiss_both_inc -> GetBinError  (i);
                error_both = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
                h1_det_eff_both -> SetBinError(i,error_both);
        }
	
	h1_det_eff_eppn     -> SetTitle("^{3}He(e,e'ppn)/^{3}He(e,e'pp)n"                             );
	h1_det_eff_eppnpipi -> SetTitle("^{3}He(e,e'pp#pi^{+}#pi^{-}n)/^{3}He(e,e'pp#pi^{+}#pi^{-})n" );
	h1_det_eff_both     -> SetTitle("Combined neutron detection efficiency"                       );

	h1_det_eff_eppn     -> GetYaxis() -> SetTitle("neut. det. eff.");
        h1_det_eff_eppnpipi -> GetYaxis() -> SetTitle("neut. det. eff.");
        h1_det_eff_both     -> GetYaxis() -> SetTitle("neut. det. eff.");

	h1_det_eff_eppn     -> SetFillColor(10);
	h1_det_eff_eppnpipi -> SetFillColor(10);
	h1_det_eff_both     -> SetFillColor(10);

	// --------------------------------------------------------------------------
	// Drawing histograms
	TCanvas * c1 = new TCanvas("c1","3He(e,e'pp)",1200,900);
	c1 -> Divide(3,3);
	c1 -> cd(1);	h1_M_miss_eppn      -> Draw();	h1_M_miss_eppn_cut       -> Draw("same");
	c1 -> cd(2);	h1_max_vz_sig_eppn  -> Draw();	h1_max_vz_sig_eppn_cut   -> Draw("same");
	c1 -> cd(3);	h1_P_miss_eppn      -> Draw();	h1_P_miss_eppn_cut       -> Draw("same");
	c1 -> cd(4);	h1_u_eppn           -> Draw();	h1_u_eppn_cut            -> Draw("same");
	c1 -> cd(5);	h1_v_eppn           -> Draw();	h1_v_eppn_cut            -> Draw("same");
	c1 -> cd(6);	h1_w_eppn           -> Draw();	h1_w_eppn_cut            -> Draw("same");
	c1 -> cd(7);    h1_cos_miss_ec_eppn -> Draw();	h1_cos_miss_ec_eppn_cut  -> Draw("same");
	c1 -> cd(8);	h1_beta_eppn        -> Draw();	h1_beta_eppn_cut         -> Draw("same");

	TCanvas * c2 = new TCanvas("c2","3He(e,e'pp pi+pi-)",1200,900);
	c2 -> Divide(3,3);
	c2 -> cd(1);    h1_M_miss_eppnpipi      -> Draw();	h1_M_miss_eppnpipi_cut       -> Draw("same");
	c2 -> cd(2);    h1_max_vz_sig_eppnpipi  -> Draw();	h1_max_vz_sig_eppnpipi_cut   -> Draw("same");
	c2 -> cd(3);    h1_P_miss_eppnpipi      -> Draw();	h1_P_miss_eppnpipi_cut       -> Draw("same");
	c2 -> cd(4);	h1_u_eppnpipi           -> Draw();	h1_u_eppnpipi_cut            -> Draw("same");
	c2 -> cd(5);	h1_v_eppnpipi           -> Draw();	h1_v_eppnpipi_cut            -> Draw("same");
	c2 -> cd(6);	h1_w_eppnpipi           -> Draw();	h1_w_eppnpipi_cut            -> Draw("same");
	c2 -> cd(7);    h1_cos_miss_ec_eppnpipi -> Draw();	h1_cos_miss_ec_eppnpipi_cut  -> Draw("same");
	c2 -> cd(8);    h1_beta_eppnpipi        -> Draw();	h1_beta_eppnpipi_cut         -> Draw("same");

	TCanvas * c3 = new TCanvas("c3","Neutron detection efficiency",1200,900);
        c3 -> Divide(3,3);
	c3 -> cd(1);	h1_pmiss_epp        -> Draw();
	c3 -> cd(2);    h1_pmiss_eppn       -> Draw();
	c3 -> cd(3);    h1_det_eff_eppn     -> Draw("BARE");
	fit_gfunction( h1_det_eff_eppn,0.5,2.1);

        c3 -> cd(4); 	h1_pmiss_epppipi    -> Draw();
	c3 -> cd(5);    h1_pmiss_epppipin   -> Draw();
	c3 -> cd(6);    h1_det_eff_eppnpipi -> Draw("BARE");
	fit_gfunction( h1_det_eff_eppnpipi,0.5,2.1);

        c3 -> cd(7); 	h1_pmiss_both_inc   -> Draw();
        c3 -> cd(8);	h1_pmiss_both_exc   -> Draw();
	c3 -> cd(9);	h1_det_eff_both     -> Draw("BARE");
	fit_gfunction( h1_det_eff_both,0.5,2.1);


	return 0;
}
// ====================================================================================================================================================
double max_outta_three ( double A , double B, double C ){

	double temp_max, max_val;

	temp_max = TMath::Max(A,B       );
	max_val  = TMath::Max(C,temp_max);

	return max_val;
}

// ====================================================================================================================================================
Double_t func_fit(Double_t *x, Double_t *par){
        Double_t xx = x[0];
        Double_t value;
        if     (xx <  par[0]){value = par[1]*xx    +par[2];}
        else if(xx >= par[0]){value = par[1]*par[0]+par[2];}
	return value;
}
// ====================================================================================================================================================
void fit_gfunction(TH1F * gPlot, float minx, float maxx){
        Double_t par[3] = {1.5,0.3,0.1};
        TF1* Func = new TF1("Func", func_fit, 0.5,2.1,3);
        Func -> SetParameters(par);
	Func -> SetParLimits(0,1.3,1.7);
        gPlot -> Fit("Func","EM","", minx, maxx);
}


