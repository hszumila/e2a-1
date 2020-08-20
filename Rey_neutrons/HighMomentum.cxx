// ==================================
// Code for the 4He/3He/12C analysis
// QE High-Momentum nucleons
// by Reynier Cruz Torres
// ==================================
#include "Riostream.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TRint.h"
#include "TLatex.h"
#include "TChain.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TRandom.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
using namespace std;

const double Me      = .000511;         	// Electron mass in GeV
const double Mp      = 0.93827;         	// Proton   mass in GeV
const double Mn      = 0.93957;         	// Neutron  mass in GeV
const double BE_3He  = 0.00772;                 // 3He binding energy in GeV
const double BE_4He  = 0.02830;                 // 4He binding energy in GeV
const double BE_12C  = 0.09215;                 // 12C binding energy in GeV
const double M3He    = 2* Mp+Mn -BE_3He;        // 3He mass in GeV
const double M4He    = 2*(Mp+Mn)-BE_4He;        // 4He mass in GeV
const double M12C    = 6*(Mp+Mn)-BE_12C;        // 12C mass in GeV
const double pi      = 3.14159;
// Hadron Correction Data
double  fgPar_Pfidft1l[6][6]   , fgPar_Pfidft1r[6][6], fgPar_Pfidft2l[6][6], fgPar_Pfidft2r[6][6], fgPar_Pfidbt1l[6][6],
	fgPar_Pfidbt1r[6][6]   , fgPar_Pfidbt2l[6][6], fgPar_Pfidbt2r[6][6], fgPar_Pfidbl  [6][6], fgPar_Pfidbr  [6][6],	
	fgPar_Pfid_For[6][4][7], fgPar_Pfid_Bak[6][4][7], fgPar_Pfid_ScpdS2[2][6],
	fgPar_Pfid_ScpdS3[8][6], fgPar_Pfid_ScpdS4[4][6], fgPar_Pfid_ScpdS5[8][6];

// Neutron momentum resolution parameters
double ec_time_res, ec_extr_res, ec_extr_res_err;

// Neutron detection efficiency parameters
TGraph *n_det_eff, *n_det_eff_err;

// Proton/Neutrin cross section parameters
TF2 *sigma_p_GK05, *sigma_p_AMT, *sigma_p_MMD, *sigma_n_GK05, *sigma_n_AMT, *sigma_n_MMD;

// Forward-declaring functions
void prep_pink_band( double GK05_res , double AMT_res , double MMD_res , double * y_ctr_array , double * y_err_array );
void load_cross_section_maps();								// Function to load proton/neutron cross sections
void pretty_TH2F(TH2F * h2);                                                            // Makes TH2F's look nice
void pretty_TH1F(TH1F * h1);                                                            // Makes TH1F's look nice
void pretty_TGraphErrors(TGraphErrors * gP);                                            // Makes TGraphError's look nice
float Coulomb_Corr_MeV(int Z, int A);							// Calculates proton Coulomb Energy Correction
TLorentzVector P4_Coulomb_Corrected(TLorentzVector P4_Uncorrected, int Z, int A);	// Applies proton Coulomb Energy Correction to proton momentum
void read_n_mom_res_params();								// Reads in parameters needed for neutron momentum resolution
void read_n_det_eff_params();								// Reads in parameters needed for neutron detection effocoency
double n_efficiency     ( double mom_GeV );                                             // Returns neutron detection efficiency
double n_efficiency_err ( double mom_GeV );                                             // Returns neutron detection efficiency error
double f_delta_p(double p);								// Returns neutron momentum resolution in units of GeV
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);	// Returns missing energy in units of GeV
bool p_points_to_EC_fid(TVector3 pm, double vtx_z, double dist);			// Checks whether a given momentum points to the EC or not
bool read_p_fid_params();								// Reads parameters for proton/neutron fiducial cuts
bool pFiducialCut(TVector3 momentum);							// Checks whether a proton/neutron is in the fiducial region

// ====================================================================================================================================================================
//Main Function
int main(int argc, char ** argv) {

	if (argc != 2){
		cerr << "\n\tWrong number of arguments. Instead try:\n\n"
			<< "\t\tneutron_ana /path/to/input/file \n\n"
			<< "\twhere the input file is a .dat file containing\n"
			<< "\ta list of addresses to all the input root files\n"
			<< "\tthat are to be analyzed\n\n"
			<< "WARNING: 1) Make sure you don't include the '.dat' extension of the filename!\n"
			<< "         2) Make sure all the root files you list in the input file correspond\n"
			<< "            to the same target and beam energy!\n\n";
		return -1;
	}

	TString filename = argv[1];
	filename.Append(".dat");	
	// ----------------------------------------------------------

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	// ----------------------------------------------------------
	// Prepare the list of root files to be imported
	int n_lines = 0;
	TString input_root[1000];
	ifstream inputfile;
	inputfile.open(filename);
	cout << "***** Will be working with the following root files *****\n";
	while(!(inputfile.eof())){
		inputfile >> input_root[n_lines];
		cout << input_root[n_lines] << endl;
		n_lines++;
	}
	n_lines--;
	inputfile.close();
	cout << "*********************************************************\n";
	// ----------------------------------------------------------
	// Load all necessary parameters
	cout << "Loading all necessary parameters" << endl;
	read_p_fid_params      ();
	read_n_det_eff_params  ();
	read_n_mom_res_params  ();
	load_cross_section_maps();
	cout << "*********************************************************\n";
	// ----------------------------------------------------------
	const double ec_theta = 65.*pi/180.; // degrees
	const double z0       = 563.1; // cm
	// ----------------------------------------------------------
	//Define Branch Variables
	int nRun, nParticles, nProtons, nNeutrons, Part_type[5];
	double Nu, Q2, mom_x[5], mom_y[5], mom_z[5], vtx_z[5], ec_x[5], ec_y[5];
	// -----------------------
	//Read in data files
	double Ebeam = 4.461; // GeV
	TChain *T = new TChain("T");	
	for(int i = 0 ; i < n_lines ; i++) T->Add(input_root[i]);
	// -----------------------
	//Define Branch Addresses
	T->SetBranchAddress("nRun"      ,&nRun      );
	T->SetBranchAddress("nParticles",&nParticles);
	T->SetBranchAddress("nProtons"  ,&nProtons  );
	T->SetBranchAddress("nNeutrons" ,&nNeutrons );
	T->SetBranchAddress("Nu"        ,&Nu        );
	T->SetBranchAddress("Q2"        ,&Q2        );
	T->SetBranchAddress("mom_x"     , mom_x     );
	T->SetBranchAddress("mom_y"     , mom_y     );
	T->SetBranchAddress("mom_z"     , mom_z     );
	T->SetBranchAddress("Part_type" , Part_type );
	T->SetBranchAddress("vtx_z_cor" , vtx_z     );
	T->SetBranchAddress("ec_x"      , ec_x      );
	T->SetBranchAddress("ec_y"      , ec_y      );
	// --------------------------------------------------------------------------------------------------
	// Obtaining run number and other important parameters
	T->GetEvent(0);
	int tab_run, tab_E1, tab_torus, tab_mini;
	TString tab_targ;
	ifstream run_table;
	run_table.open("./parameters/run_table.dat");
	do{
		run_table >> tab_run  ;
		run_table >> tab_E1   ;
		run_table >> tab_torus;
		run_table >> tab_mini ;
		run_table >> tab_targ ;
	}
	while(tab_run != nRun);

	cout<<"Ebeam  = "<< tab_E1<<"\nTorus  = "<<tab_torus<<"\nMini   = "<<tab_mini<<"\nTarget = "<<tab_targ<<endl;
	if      (tab_E1 != 4461){cout << "Error: This code was created for 4.4 GeV data" << endl;	exit(0);}

	// ----------------------------------------------------------
	//Define Variables
	int nEnt_ep;

	double  Pmiss, Emiss, Mmiss, p_ecx, p_ecy, p_smeared, Pmiss_sm, Emiss_sm, Mmiss_sm, 
		Xb, r_ec, pm_theta, third_theta, theta, phi, theta_e,
		Z, N, theta_pq, p_over_q, psm_over_q, n_weight;
	// --------------------
	// Target Mass
	double Mtar = -999.;
	if     (tab_targ == "4He"){	Mtar = M4He;	Z = 2;	N = 2;	}
	else if(tab_targ == "3He"){	Mtar = M3He;	Z = 2;	N = 1;	}
	else if(tab_targ == "12C"){	Mtar = M12C;	Z = 6;	N = 6;	}
	// --------------------
	TVector3 u1;
	TLorentzVector V4_Targ (0,0,0    ,Mtar );
	TLorentzVector V4_E0   (0,0,Ebeam,Ebeam);
	TLorentzVector q4, e_mom, h_mom, p_sm_mom;	
	TLorentzVector p_pair  (0,0,0,2*Mp);
	TLorentzVector n_pair  (0,0,0,2*Mn);
	// --------------------------------------------------
	double Pmiss_low  = 0.30; // GeV
	double Pmiss_high = 1.00; // GeV
	double Mmiss_high = 1.10; // GeV
	double th_Nq_low  =  0.0; // deg
	double th_Nq_high = 25.0; // deg
	double p_p_q_low  = 0.62; // GeV
	double p_p_q_high = 0.96; // GeV
	double n_p_q_low  = 0.62; // GeV
	double n_p_q_high = 1.10; // GeV
	// --------------------------------------------------
	// Optimization of QE cuts
	const int no_mmiss_cuts = 7;
	float Mmiss_cut_lim[no_mmiss_cuts] = {0};
	for(int i = 0 ; i < no_mmiss_cuts ; i++) Mmiss_cut_lim[i] = 1.101+0.037*(float)i;

	const int no_pmiss_cuts = 25;
	float Pmiss_cut_lim[no_pmiss_cuts] = {0};
	for(int i = 0 ; i < no_pmiss_cuts ; i++) Pmiss_cut_lim[i] = 0.01+0.0245*(float)i;

	float fal_denom  [no_mmiss_cuts][no_pmiss_cuts] = {{0}};

	float fal_pos_num[no_mmiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_pos    [no_mmiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_pos_e  [no_mmiss_cuts][no_pmiss_cuts] = {{0}};

	float fal_neg_num[no_mmiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_neg    [no_mmiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_neg_e  [no_mmiss_cuts][no_pmiss_cuts] = {{0}}; 

	// --------------------
	// Variables for cross section ratios
	double nu_eep = 0;
	double nu_een = 0;
	double raw_nu_p = 0;
	double raw_nu_n = 0;
	// --------------------
        // Variables for theoretical cross section ratios
        double nu_eep_een_GK05 = 0;
        double nu_eep_een_AMT  = 0;
        double nu_eep_een_MMD  = 0; 
	// ----------------------------------------------------------
	//Define Histograms

	// Acceptance matching: 10 cm EC cut
	TH2F * h2_en_EC_xy       = new TH2F("h2_en_EC_xy"      ,";EC_{x} [cm];EC_{y} [cm]"      , 100 , -500, 500 , 100 , -500, 500 );
	TH2F * h2_ep_EC_xy       = new TH2F("h2_ep_EC_xy"      ,";EC_{x} [cm];EC_{y} [cm]"      , 100 , -500, 500 , 100 , -500, 500 );
	TH2F * h2_ep_EC_xy_cut   = new TH2F("h2_ep_EC_xy_cut"  ,";EC_{x} [cm];EC_{y} [cm]"      , 100 , -500, 500 , 100 , -500, 500 );

	TH1F * h1_en_EC_x        = new TH1F("h1_en_EC_x"       ,"EC x distribution;EC_{x} [cm];Counts" , 100 , -500, 500 );
	TH1F * h1_ep_EC_x        = new TH1F("h1_ep_EC_x"       ,"EC x distribution;EC_{x} [cm];Counts" , 100 , -500, 500 );
	TH1F * h1_ep_EC_x_cut    = new TH1F("h1_ep_EC_x_cut"   ,"EC x distribution;EC_{x} [cm];Counts" , 100 , -500, 500 );

	TH1F * h1_en_EC_y        = new TH1F("h1_en_EC_y"       ,"EC y distribution;EC_{y} [cm];Counts" , 100 , -500, 500 );
	TH1F * h1_ep_EC_y        = new TH1F("h1_ep_EC_y"       ,"EC y distribution;EC_{y} [cm];Counts" , 100 , -500, 500 );
	TH1F * h1_ep_EC_y_cut    = new TH1F("h1_ep_EC_y_cut"   ,"EC y distribution;EC_{y} [cm];Counts" , 100 , -500, 500 );

	h2_ep_EC_xy_cut -> SetMarkerColor(2);
	h2_ep_EC_xy_cut -> SetLineColor  (2);

	h1_ep_EC_x	-> SetLineColor(2);
	h1_ep_EC_x_cut  -> SetLineColor(8);
	h1_ep_EC_y      -> SetLineColor(2);
	h1_ep_EC_y_cut  -> SetLineColor(8);
	// ----------------------------------
	// Acceptance matching: phi-theta cut
	TH2F * h2_ep_th_phi      = new TH2F("h2_ep_th_phi"      ,";#phi [deg];#theta [deg]"    ,100,-100.,390.,100,0.,70.);
	TH2F * h2_en_th_phi      = new TH2F("h2_en_th_phi"      ,";#phi [deg];#theta [deg]"    ,100,-100.,390.,100,0.,70.);
	TH2F * h2_en_th_phi_cut  = new TH2F("h2_en_th_phi_cut"  ,";#phi [deg];#theta [deg]"    ,100,-100.,390.,100,0.,70.);

	TH1F * h1_en_th          = new TH1F("h1_en_th"          ,";#theta [deg];Counts"               ,100,   0.,55. );
	TH1F * h1_ep_th          = new TH1F("h1_ep_th"          ,";#theta [deg];Counts"               ,100,   0.,55. );
	TH1F * h1_en_th_cut      = new TH1F("h1_en_th_cut"      ,";#theta [deg];Counts"               ,100,   0.,55. );

	TH1F * h1_en_phi         = new TH1F("h1_en_phi"         ,";#phi [deg];Counts"                 ,100,-150.,390.);
	TH1F * h1_ep_phi         = new TH1F("h1_ep_phi"         ,";#phi [deg];Counts"                 ,100,-150.,390.);
	TH1F * h1_en_phi_cut     = new TH1F("h1_en_phi_cut"     ,";#phi [deg];Counts"                 ,100,-150.,390.);

	h2_en_th_phi_cut -> SetMarkerColor(2);
	h2_en_th_phi_cut -> SetLineColor  (2);

	h1_en_th       -> SetLineColor(2);
	h1_en_th_cut   -> SetLineColor(8);
	h1_en_phi      -> SetLineColor(2);
	h1_en_phi_cut  -> SetLineColor(8);
	// ----------------------------------
	// Momentum cut
	TH1F * h1_n_p            = new TH1F("h1_n_p"             ,";p [GeV];Counts"                   ,100,   0.,2.8 );
	TH1F * h1_n_p_cut        = new TH1F("h1_n_p_cut"         ,";p [GeV];Counts"                   ,100,   0.,2.8 );
	TH1F * h1_p_p            = new TH1F("h1_p_p"             ,";p [GeV];Counts"                   ,100,   0.,2.8 ); 
	TH1F * h1_p_p_cut        = new TH1F("h1_p_p_cut"         ,";p [GeV];Counts"                   ,100,   0.,2.3 );
	h1_p_p     -> SetLineColor( 2);
	h1_p_p_cut -> SetLineColor( 8);
	h1_n_p_cut -> SetLineColor(93);
	// ----------------------------------
	TH1F * h1_Q2_Xb_cut_p   = new TH1F("h1_Q2_Xb_cut_p"      ,"protons (x_{B}>1.2 cut);Q^{2} [GeV^{2}];Counts"                      , 25 , 0.5 , 5.6 );
	TH1F * h1_Q2_Xb_cut_psm = new TH1F("h1_Q2_Xb_cut_psm"    ,"smeared protons and neutrons (x_{B}>1.1 cut);Q^{2} [GeV^{2}];Counts" , 15 , 0.5 , 5.6 );
	TH1F * h1_Q2_Xb_cut_n   = new TH1F("h1_Q2_Xb_cut_n"      ,"smeared protons and neutrons (x_{B}>1.1 cut);Q^{2} [GeV^{2}];Counts" , 15 , 0.5 , 5.6 );
	h1_Q2_Xb_cut_psm -> SetLineColor(2);

	TH2F * h2_thpq_ppq   = new TH2F("h2_thpq_ppq"  ,"protons (acc + x_{B} + p_{miss} cuts);p_{p}/q;#theta_{pq}",60,0.1,1.5,50,0,80);
	TH2F * h2_thpq_npq   = new TH2F("h2_thpq_npq"  ,"neutrons (acc + x_{B} cuts);p_{n}/q;#theta_{nq}"          ,40,0.1,1.5,40,0,80);
	TH2F * h2_thpq_psmpq = new TH2F("h2_thpq_psmpq","smeared protons (acc + x_{B} cuts);p_{p}/q;#theta_{pq}"   ,40,0.1,1.5,40,0,80);

	TH2F * h2_thnq_pn    = new TH2F("h2_thpq_pn"   ,"neutrons (acc + x_{B} cuts);#theta_{pq};p_{n} [GeV]"    ,40,0,80,40,0.7,2.5);
	TH2F * h2_thnq_thn   = new TH2F("h2_thpq_thn"  ,"neutrons (acc + x_{B} cuts);#theta_{pq};#theta_{n}"     ,40,0,80,40,18.,45.);
	TH2F * h2_thnq_pmiss = new TH2F("h2_thpq_pmiss","neutrons (acc + x_{B} cuts);#theta_{pq};p_{miss} [GeV]" ,40,0,80,40,0.0,1.4);
	TH2F * h2_thnq_Q2    = new TH2F("h2_thpq_Q2"   ,"neutrons (acc + x_{B} cuts);#theta_{pq};Q^{2} [GeV^{2}]",40,0,80,40,1.0,5.0);
	TH2F * h2_thnq_Nu    = new TH2F("h2_thpq_Nu"   ,"neutrons (acc + x_{B} cuts);#theta_{pq};#omega [GeV]"   ,40,0,80,40,0.4,2.3);
	TH2F * h2_thnq_Xb    = new TH2F("h2_thpq_Xb"   ,"neutrons (acc + x_{B} cuts);#theta_{pq};x_{B}"          ,40,0,80,40,1.0,2.2);

	TH1F * h1_Mmiss_p    = new TH1F("h1_Mmiss_p"   ,"protons (acc + x_{B} + p_{miss} + leading hadron cuts);m_{miss} [GeV];Counts"           ,40,0,2);
	TH1F * h1_Mmiss_psm  = new TH1F("h1_Mmiss_psm" ,"smeared protons and neutrons (acc + x_{B} + leading hadron cuts);m_{miss} [GeV];Counts" ,25,0,2);
	TH1F * h1_Mmiss_n    = new TH1F("h1_Mmiss_n"   ,"smeared protons and neutrons (acc + x_{B} + leading hadron cuts);m_{miss} [GeV];Counts" ,25,0,2);
	h1_Mmiss_psm -> SetLineColor(2);

	TH1F * h1_Pmiss_p    = new TH1F("h1_Pmiss_p"   ,"protons (acc + x_{B} + p,m_{miss} + leading hadron cuts);p_{miss} [GeV];Counts"          ,40,0,1.4);
	TH1F * h1_Pmiss_psm  = new TH1F("h1_Pmiss_psm" ,"smeared protons and neutrons (acc + x_{B} + leading hadron cuts);p_{miss} [GeV];Counts"  ,20,0,1.4);
	TH1F * h1_Pmiss_n    = new TH1F("h1_Pmiss_n"   ,"smeared protons and neutrons (acc + x_{B} + leading hadron cuts);p_{miss} [GeV];Counts"  ,20,0,1.4);
	h1_Pmiss_psm -> SetLineColor(2);

	TH1F * h1_Pm_psm     = new TH1F("h1_Pm_psm"   ,";p_{miss} [GeV];Counts",40,0,1.4);
	TH1F * h1_Pm_psm_fp  = new TH1F("h1_Pm_psm_fp",";p_{miss} [GeV];Counts",40,0,1.4);
	TH1F * h1_Pm_psm_fn  = new TH1F("h1_Pm_psm_fn",";p_{miss} [GeV];Counts",40,0,1.4);
	h1_Pm_psm    -> SetLineColor(2);
	h1_Pm_psm_fp -> SetLineColor(1);
	h1_Pm_psm_fn -> SetLineColor(8);
	// ----------------------------------
	// Missing energy distribution
	TH1F * h1_Emiss_p   = new TH1F("h1_Emiss_p"  ,"protons (acc + x_{B} + p,m_{miss} + leading hadron cuts);E_{miss} [GeV];Counts"                           ,40,-.3,.7);
	TH1F * h1_Emiss_psm = new TH1F("h1_Emiss_psm","smeared protons and neutrons (acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);E_{miss} [GeV];Counts",20,-.3,.7);
	TH1F * h1_Emiss_n   = new TH1F("h1_Emiss_n"  ,"smeared protons and neutrons (acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);E_{miss} [GeV];Counts",20,-.3,.7);
	h1_Emiss_psm -> SetLineColor(2);
	// ----------------------------------
	// Electron kinematics variables
	TH1F * h1_pe_psm = new TH1F("h1_pe_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{e} [GeV];Counts"     ,20,2.6,4.5);
	TH1F * h1_pe_n   = new TH1F("h1_pe_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{e} [GeV];Counts"     ,20,2.6,4.5);
	TH1F * h1_te_psm = new TH1F("h1_te_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{e} [deg];Counts",20,0.0,40.);
	TH1F * h1_te_n   = new TH1F("h1_te_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{e} [deg];Counts",20,0.0,40.);
	TH1F * h1_Q2_psm = new TH1F("h1_Q2_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);Q^{2} [GeV^{2}];Counts" ,20,1.0,5.0);
	TH1F * h1_Q2_n   = new TH1F("h1_Q2_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);Q^{2} [GeV^{2}];Counts" ,20,1.0,5.0);
	TH1F * h1_Nu_psm = new TH1F("h1_Nu_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#omega [GeV];Counts"    ,20,0.3,2.3);
	TH1F * h1_Nu_n   = new TH1F("h1_Nu_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#omega [GeV];Counts"    ,20,0.3,2.3);
	TH1F * h1_Xb_psm = new TH1F("h1_Xb_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);x_{B};Counts"           ,20,1.0,2.2);
	TH1F * h1_Xb_n   = new TH1F("h1_Xb_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);x_{B};Counts"           ,20,1.0,2.2);
	h1_pe_psm -> SetLineColor(2);
	h1_te_psm -> SetLineColor(2);
	h1_Q2_psm -> SetLineColor(2);
	h1_Nu_psm -> SetLineColor(2);
	h1_Xb_psm -> SetLineColor(2);
	// ----------------------------------
	// Proton, neutron kinematics variables
	TH1F * h1_pm_psm  = new TH1F("h1_pm_psm" ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{miss} [GeV];Counts"   ,20,0.3,1.2);
	TH1F * h1_pm_n    = new TH1F("h1_pm_n"   ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{miss} [GeV];Counts"   ,20,0.3,1.2);
	TH1F * h1_ph_psm  = new TH1F("h1_ph_psm" ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{N} [GeV];Counts"      ,20,0.7,2.5);
	TH1F * h1_ph_n    = new TH1F("h1_ph_n"   ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);p_{N} [GeV];Counts"      ,20,0.7,2.5);
	TH1F * h1_th_psm  = new TH1F("h1_th_psm" ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{N} [deg];Counts" ,20,18.,45.);
	TH1F * h1_th_n    = new TH1F("h1_th_n"   ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{N} [deg];Counts" ,20,18.,45.);
	TH1F * h1_thq_psm = new TH1F("h1_thq_psm","(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{Nq} [deg];Counts", 5, 0.,25.);
	TH1F * h1_thq_n   = new TH1F("h1_thq_n"  ,"(acc + x_{B} + p,m^{opt}_{miss} + leading hadron cuts);#theta_{Nq} [deg];Counts", 5, 0.,25.);
	h1_pm_psm  -> SetLineColor(2);
	h1_ph_psm  -> SetLineColor(2);
	h1_th_psm  -> SetLineColor(2);
	h1_thq_psm -> SetLineColor(2);	

	// -------------------------------------------------------------------------------------------------------------------------------------
	//Process Data
	nEnt_ep = T->GetEntries();	
	// --------------------------------------------------------------------------------------
	cout << "**************************\nProcessing (e,e'p) events\n";
	for(int evt=0;evt<nEnt_ep;evt++){
		if(evt%1000000==0) cout << " events processed = " << evt << " out of " << nEnt_ep << endl;
		T->GetEntry(evt);	
		if(nProtons!=1) continue;	
		if((Part_type[0]!=-11)||(Part_type[1]!=2212)){cout << "Error: event in (e,e'p) section doesn't have correct content of particles" << endl; exit(1);}

		e_mom.SetXYZM(mom_x[0],mom_y[0],mom_z[0],Me);	// recoiling electron momentum
		h_mom.SetXYZM(mom_x[1],mom_y[1],mom_z[1],Mp);	// struck proton momentum
		h_mom = P4_Coulomb_Corrected(h_mom,Z,Z+N);      // Applying Coulomb corrections

		q4    = V4_E0-e_mom;

		Pmiss = (h_mom - q4).Rho();			// Missing momentum
		Emiss = fn_Emiss(Pmiss,Nu,Mtar,h_mom.E(),Mp);	// Missing energy

		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		// EC 10 cm cut
		pm_theta    = h_mom.Vect().Theta();
		third_theta = (pi - ec_theta - pm_theta);
		r_ec = (z0-vtx_z[1])*sin(ec_theta)/sin(third_theta);
		p_ecx = r_ec*sin(h_mom.Vect().Theta())*cos(h_mom.Vect().Phi());
		p_ecy = r_ec*sin(h_mom.Vect().Theta())*sin(h_mom.Vect().Phi());
		h1_ep_EC_x  -> Fill(p_ecx       );
		h1_ep_EC_y  -> Fill(p_ecy       );
		h2_ep_EC_xy -> Fill(p_ecx,p_ecy );
		if(!p_points_to_EC_fid(h_mom.Vect(),vtx_z[0],10.)) continue;
		h1_ep_EC_x_cut    -> Fill(p_ecx      );
		h1_ep_EC_y_cut    -> Fill(p_ecy      );
		h2_ep_EC_xy_cut   -> Fill(p_ecx,p_ecy);
		// ---------------------------------------
		// phi-theta fiducial cut
		theta = h_mom.Vect().Theta()*180./pi;
		phi   = h_mom.Vect().Phi()  *180./pi;	if(phi<-30) phi+=360;

		if(theta>44.) continue;	// Cutting top of EC
		h1_ep_th     -> Fill(    theta);
		h1_ep_phi    -> Fill(      phi);
		h2_ep_th_phi -> Fill(phi,theta);
		// ---------------------------------------
		// Momentum cut
		h1_p_p     -> Fill(h_mom.Vect().Mag());
		if(h_mom.Vect().Mag()>2.34) continue;
		h1_p_p_cut -> Fill(h_mom.Vect().Mag());
		// ------------------------------------------------------------------------------ 

		Mmiss    = (q4+p_pair-h_mom).Mag();

		// ------------------------------------------------------------------------------ 
		// Smearing protons to simulate neutrons
		p_smeared = gRandom -> Gaus(h_mom.Vect().Mag(),f_delta_p(h_mom.Vect().Mag()));
		if(p_smeared>2.34) continue;

		// Copying direction of p momentum vector
		u1 = h_mom.Vect().Unit();
		p_sm_mom.SetXYZM(p_smeared*u1.X(),p_smeared*u1.Y(),p_smeared*u1.Z(),Mp);

		Pmiss_sm = (p_sm_mom - q4).Rho();		// Missing momentum
		Mmiss_sm = (q4+p_pair-p_sm_mom).Mag(); 		// Missing mass

		Xb = Q2/(2.*Nu*Mp);
		if(Xb<1.1) continue;
		if(h_mom   .Rho()<0.5) continue;
		if(p_sm_mom.Rho()<0.5) continue;

		u1 = h_mom.Vect().Unit();
		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;
		p_over_q = h_mom.Rho()/q4.Rho();
		psm_over_q = p_sm_mom.Rho()/q4.Rho();
		// ------------------------------------------------------------------------------------------
		// False positive and False negative rates
		if((theta_pq < th_Nq_high )&&(theta_pq > th_Nq_low )&&(psm_over_q > n_p_q_low)&&(psm_over_q < n_p_q_high)){
			for(int i = 0 ; i < no_pmiss_cuts ; i++){
				for(int j = 0 ; j < no_mmiss_cuts ; j++){
					// -------
					// False positive
					if((Pmiss_sm>Pmiss_cut_lim[i])&&(Pmiss_sm<Pmiss_high)&&(Mmiss_sm<Mmiss_cut_lim[j])){ 

						fal_denom[j][i]++;
						if(!((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)&&(Mmiss<Mmiss_high))){
							fal_pos_num[j][i]++;
						}

					}
					// -------
					// False negative
					if((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)&&(Mmiss<Mmiss_high)){
						if(!((Pmiss_sm>Pmiss_cut_lim[i])&&(Pmiss_sm<Pmiss_high)&&(Mmiss_sm<Mmiss_cut_lim[j]))){
							fal_neg_num[j][i]++;
						}
					}
				}
			}
		}
	}

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Doing other calculations
	// False positive rate
	double A, B, eA, eB;

	TGraphErrors ** f_fal_pos = new TGraphErrors*[no_mmiss_cuts];
	TLegend * leg_false_pos = new TLegend(0.65,0.50,0.85,0.85);
	leg_false_pos -> SetLineColor(0);
	leg_false_pos -> SetHeader("M_{miss} cut [GeV]");
	for(int j = 0 ; j < no_mmiss_cuts ; j++){
		for(int i = 0 ; i < no_pmiss_cuts ; i++){
			A  = (float)fal_pos_num[j][i];
			B  = (float)fal_denom  [j][i];
			eA = sqrt(A);
			eB = sqrt(B);
			fal_pos  [j][i] = 100.*A/B;
			fal_pos_e[j][i] = 0.;//100.*sqrt(eA*eA - A*A*eB*eB/(B*B))/B;
		}
		f_fal_pos[j] = new TGraphErrors(no_pmiss_cuts,Pmiss_cut_lim,fal_pos[j],0,fal_pos_e[j]);
		f_fal_pos[j] -> SetMarkerColor(j+1);
		f_fal_pos[j] -> SetMarkerSize (1  );
		f_fal_pos[j] -> SetMarkerStyle(20 );
		f_fal_pos[j] -> SetFillColor  (0  );
		f_fal_pos[j] -> SetMinimum(  0);
		f_fal_pos[j] -> SetMaximum( 60);
		f_fal_pos[j] -> SetTitle("");
		f_fal_pos[j] -> GetXaxis() -> SetTitle("p_{miss} cut [GeV]");
		f_fal_pos[j] -> GetYaxis() -> SetTitle("False Positive [%]");
		leg_false_pos -> AddEntry(f_fal_pos[j],Form("%.3f",Mmiss_cut_lim[j]));
	}
	// False negative rate
	TGraphErrors ** f_fal_neg = new TGraphErrors*[no_mmiss_cuts];
	TLegend * leg_false_neg = new TLegend(0.20,0.50,0.40,0.85);
	leg_false_neg -> SetLineColor(0);
	leg_false_neg -> SetHeader("M_{miss} cut [GeV]");
	for(int j = 0 ; j < no_mmiss_cuts ; j++){
		for(int i = 0 ; i < no_pmiss_cuts ; i++){
			A  = (float)fal_neg_num[j][i];
			B  = (float)fal_denom  [j][i];
			eA = sqrt(A);
			eB = sqrt(B);
			fal_neg  [j][i] = 100.*A/B;
			fal_neg_e[j][i] = 0.;//100.*sqrt(eA*eA - A*A*eB*eB/(B*B))/B;
		}
		f_fal_neg[j] = new TGraphErrors(no_pmiss_cuts,Pmiss_cut_lim,fal_neg[j],0,fal_neg_e[j]);
		f_fal_neg[j] -> SetMarkerColor(j+1);
		f_fal_neg[j] -> SetMarkerSize (1  );
		f_fal_neg[j] -> SetMarkerStyle(20 );
		f_fal_neg[j] -> SetFillColor  (0  );
		f_fal_neg[j] -> SetMinimum(  0);
		f_fal_neg[j] -> SetMaximum(100);
		f_fal_neg[j] -> SetTitle("");
		f_fal_neg[j] -> GetXaxis() -> SetTitle("p_{miss} cut [GeV]");
		f_fal_neg[j] -> GetYaxis() -> SetTitle("False Negative [%]");
		leg_false_neg -> AddEntry(f_fal_neg[j],Form("%.3f",Mmiss_cut_lim[j]));
	}
	// False negative - false positive
	TGraphErrors ** f_pos_neg = new TGraphErrors*[no_mmiss_cuts];
	TLegend * leg_pos_neg = new TLegend(0.65,0.50,0.85,0.85);
	leg_pos_neg -> SetLineColor(0);
	leg_pos_neg -> SetHeader("M_{miss} cut [GeV]");
	for(int j = 0 ; j < no_mmiss_cuts ; j++){
		f_pos_neg[j] = new TGraphErrors(no_pmiss_cuts,fal_pos[j],fal_neg[j],fal_pos_e[j],fal_neg_e[j]);
		f_pos_neg[j] -> SetMarkerColor(j+1);
		f_pos_neg[j] -> SetMarkerSize (1  );
		f_pos_neg[j] -> SetMarkerStyle(20 );
		f_pos_neg[j] -> SetFillColor  (0  );
		f_pos_neg[j] -> GetXaxis() -> SetLimits(0,100);
		f_pos_neg[j] -> SetMinimum(  0);
		f_pos_neg[j] -> SetMaximum(100);
		f_pos_neg[j] -> SetTitle("");
		f_pos_neg[j] -> GetXaxis() -> SetTitle("False Positive [%]");
		f_pos_neg[j] -> GetYaxis() -> SetTitle("False Negative [%]");
		leg_pos_neg -> AddEntry(f_pos_neg[j],Form("%.3f",Mmiss_cut_lim[j]));
	}
	// Optimizing cut
	int min_i     = 99999;
	int min_j     = 99999;
	float min_val = 99999;
	float tem_val = 99999;
	for(int i = 0 ; i < no_pmiss_cuts ; i++){
		for(int j = 0 ; j < no_mmiss_cuts ; j++){
			tem_val = fal_pos[j][i]*fal_pos[j][i] + fal_neg[j][i]*fal_neg[j][i];
			if (tem_val < min_val){
				min_val = tem_val;
				min_i = i;
				min_j = j;
			}
		}
	}
	// ----
	// Forcing the values
	min_i = 16;
	min_j =  2;
	//min_i = 16;
	//min_j =  3;

	cout << "**************************\nThe simultaneous minimization of false positive and negative rate happens for:\n" <<
		"Pmiss_cut = " << Pmiss_cut_lim[min_i] << " GeV, and\nMmiss_cut = " << Mmiss_cut_lim[min_j] << " GeV\n" <<
		"and corresponds to a false positive of: " << fal_pos[min_j][min_i] << "%\nand a false negative of: "<<
		fal_neg[min_j][min_i] << "%\n";

	TLine * l_f_pos_neg = new TLine (0,0,100,100);
	l_f_pos_neg -> SetLineColor(2);

	// -------------------------------------------------------------------------------------------------------------------------------------
	cout << "**************************\nProcessing (e,e'p) events\n(a second time after QE cuts optimization)\n";
	for(int evt=0;evt<nEnt_ep;evt++){
		if(evt%1000000==0) cout << " events processed = " << evt << " out of " << nEnt_ep << endl;
		T->GetEntry(evt);
		if(nProtons!=1) continue;
		if((Part_type[0]!=-11)||(Part_type[1]!=2212)){cout << "Error: event in (e,e'p) section doesn't have correct content of particles" << endl; exit(1);}

		e_mom.SetXYZM(mom_x[0],mom_y[0],mom_z[0],Me);   // recoiling electron momentum
		h_mom.SetXYZM(mom_x[1],mom_y[1],mom_z[1],Mp);   // struck proton momentum
		h_mom = P4_Coulomb_Corrected(h_mom,Z,Z+N);      // Applying Coulomb corrections

		q4    = V4_E0-e_mom;
		Pmiss = (h_mom - q4).Rho();                     // Missing momentum
		Emiss = fn_Emiss(Pmiss,Nu,Mtar,h_mom.E(),Mp);   // Missing energy

		theta = h_mom.Vect().Theta()*180./pi;
		phi   = h_mom.Vect().Phi()  *180./pi;	if(phi  <-30) phi  +=360;		

		theta_e = e_mom.Vect().Theta()*180./pi;

		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		if(!p_points_to_EC_fid(h_mom.Vect(),vtx_z[0],10.)) continue;	// EC 10 cm cut 
		if(theta>44.) continue;						// Cutting top of EC
		// ---------------------------------------
		// Momentum cut
		if(h_mom.Vect().Mag()>2.34) continue;
		// ------------------------------------------------------------------------------	

		Xb = Q2/(2.*Nu*Mp);

		u1 = h_mom.Vect().Unit();
		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;
		p_over_q = h_mom.Rho()/q4.Rho();
		Mmiss    = (q4+p_pair-h_mom).Mag();

		if(Xb<1.2) continue;
		if(h_mom.Rho()<0.5) continue;

		h1_Q2_Xb_cut_p -> Fill( Q2 );

		if((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)){
			h2_thpq_ppq    -> Fill(p_over_q,theta_pq);
		}

		if((theta_pq < th_Nq_high )&&(theta_pq > th_Nq_low )&&(p_over_q > p_p_q_low)&&(p_over_q < p_p_q_high)){
			if((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)){
				h1_Mmiss_p -> Fill(Mmiss);

				if(Mmiss<Mmiss_high){
					h1_Pmiss_p -> Fill(Pmiss);
					h1_Emiss_p -> Fill(Emiss);

					//if(p_points_to_EC_fid((q4+V4_Targ).Vect(),vtx_z[0],10.)){
						raw_nu_p += 1.;
						nu_eep   += 1.;

						nu_eep_een_GK05 += (sigma_p_GK05->Eval(Q2,theta_e))/(sigma_n_GK05->Eval(Q2,theta_e));
	                                	nu_eep_een_AMT  += (sigma_p_AMT ->Eval(Q2,theta_e))/(sigma_n_AMT ->Eval(Q2,theta_e));
        	                        	nu_eep_een_MMD  += (sigma_p_MMD ->Eval(Q2,theta_e))/(sigma_n_MMD ->Eval(Q2,theta_e));
					//}
				}
			}
		}
	}
	// -------------------------------------------------------------------------------------------------------------------------------------
	cout << "**************************\nProcessing (e,e'p_smeared) events\n";
	for(int evt=0;evt<nEnt_ep;evt++){
		if(evt%1000000==0) cout << " events processed = " << evt << " out of " << nEnt_ep << endl;
		T->GetEntry(evt);
		if(nProtons!=1) continue;
		if((Part_type[0]!=-11)||(Part_type[1]!=2212)){cout << "Error: event in (e,e'p) section doesn't have correct content of particles" << endl; exit(1);}

		e_mom.SetXYZM(mom_x[0],mom_y[0],mom_z[0],Me);   // recoiling electron momentum
		h_mom.SetXYZM(mom_x[1],mom_y[1],mom_z[1],Mp);   // struck proton momentum
		h_mom = P4_Coulomb_Corrected(h_mom,Z,Z+N);      // Applying Coulomb corrections

		q4    = V4_E0-e_mom;
		Pmiss = (h_mom - q4).Rho();                     // Missing momentum
		Emiss = fn_Emiss(Pmiss,Nu,Mtar,h_mom.E(),Mp);   // Missing energy
		Mmiss = (q4+p_pair-h_mom).Mag();

		// ------------------------------------------------------------------------------ 
		// Smearing protons to simulate neutrons
		p_smeared = gRandom -> Gaus(h_mom.Vect().Mag(),f_delta_p(h_mom.Vect().Mag()));	
		if(p_smeared>2.34) continue;

		// Copying direction of p momentum vector
		u1 = h_mom.Vect().Unit();
		p_sm_mom.SetXYZM(p_smeared*u1.X(),p_smeared*u1.Y(),p_smeared*u1.Z(),Mp);

		Pmiss_sm = (p_sm_mom - q4).Rho();                       // Missing momentum
		Emiss_sm = fn_Emiss(Pmiss_sm,Nu,Mtar,p_sm_mom.E(),Mp);  // Missing energy
		// -----------------------
		theta = p_sm_mom.Vect().Theta()*180./pi;
		phi   = p_sm_mom.Vect().Phi()  *180./pi;	if(phi  <-30) phi  +=360;		

		theta_e = e_mom.Vect().Theta()*180./pi;
		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		if(!p_points_to_EC_fid(p_sm_mom.Vect(),vtx_z[0],10.)) continue;	// EC 10 cm cut
		if(theta>44.) continue;						// Cutting top of EC
		// ---------------------------------------
		// Momentum cut
		if(p_sm_mom.Rho()>2.34) continue;
		if(p_sm_mom.Rho()<0.5 ) continue;
		// ------------------------------------------------------------------------------

		Xb = Q2/(2.*Nu*Mp);

		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;
		p_over_q = p_sm_mom.Rho()/q4.Rho();
		Mmiss_sm = (q4+p_pair-p_sm_mom).Mag();

		if(Xb<1.1) continue;
		if(h_mom.Rho()<0.5) continue;

		h1_Q2_Xb_cut_psm -> Fill( Q2 );
		h2_thpq_psmpq    -> Fill(p_over_q,theta_pq);

		if((theta_pq < th_Nq_high )&&(theta_pq >= th_Nq_low )&&(p_over_q > n_p_q_low)&&(p_over_q < n_p_q_high)){
			h1_Mmiss_psm -> Fill(Mmiss_sm);
			h1_Pmiss_psm -> Fill(Pmiss_sm);


			if((Pmiss_sm>Pmiss_cut_lim[min_i])&&(Mmiss_sm<Mmiss_cut_lim[min_j])&&(Pmiss_sm<Pmiss_high)){
				h1_Pm_psm    -> Fill(Pmiss_sm);
				h1_Emiss_psm -> Fill(Emiss_sm);

				h1_pe_psm -> Fill(e_mom.Rho());
				h1_te_psm -> Fill(theta_e    );
				h1_Q2_psm -> Fill(Q2         );
				h1_Nu_psm -> Fill(Nu         );
				h1_Xb_psm -> Fill(Xb         );

				h1_pm_psm  -> Fill(Pmiss_sm      );
				h1_ph_psm  -> Fill(p_sm_mom.Rho());
				h1_th_psm  -> Fill(theta         );
				h1_thq_psm -> Fill(theta_pq      );

				// False positive in Pmiss distribution
				if(!((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)&&(Mmiss<Mmiss_high))){
					h1_Pm_psm_fp -> Fill(Pmiss_sm);
				}
			}
			// False negative in Pmiss distribution
			if(!((Pmiss_sm>Pmiss_cut_lim[min_i])&&(Mmiss_sm<Mmiss_cut_lim[min_j])&&(Pmiss_sm<Pmiss_high))){ 
				if((Pmiss>Pmiss_low)&&(Pmiss<Pmiss_high)&&(Mmiss<Mmiss_high)){
					h1_Pm_psm_fn -> Fill(Pmiss_sm);
				}
			}
		}

	}
	// --------------------------------------------------------------------------------------
	cout << "**************************\nProcessing (e,e'n) events\n";
	for(int evt=0;evt<nEnt_ep;evt++){
		if(evt%1000000==0) cout << " events processed = " << evt << " out of " << nEnt_ep  << endl;
		T->GetEntry(evt);
		if(nNeutrons!=1) continue;
		if((Part_type[0]!=-11)||(Part_type[2]!=2112)){cout << "Error: event in (e,e'n) section doesn't have correct content of particles" << endl; exit(1);}

		e_mom.SetXYZM(mom_x[0],mom_y[0],mom_z[0],Me);   // recoiling electron momentum
		h_mom.SetXYZM(mom_x[2],mom_y[2],mom_z[2],Mn);   // struck neutron momentum

		q4    = V4_E0-e_mom;
		Pmiss = (h_mom - q4).Rho();                     // Missing momentum
		Emiss = fn_Emiss(Pmiss,Nu,Mtar,h_mom.E(),Mn);   // Missing energy

		// ---------------------------------------
                // Weight from neutron detection eff.
                if(h_mom.Rho()<0.5) continue;
		n_weight = 1./n_efficiency(h_mom.Rho());

		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		// phi-theta fiducial cut
		theta   = h_mom.Vect().Theta()*180./pi;
		phi     = h_mom.Vect().Phi()  *180./pi;	if(phi  <-30) phi  +=360;	

		theta_e = e_mom.Vect().Theta()*180./pi;

		h1_en_th     -> Fill(    theta);
		h1_en_phi    -> Fill(      phi);
		h2_en_th_phi -> Fill(phi,theta);

		if(theta>44.) continue;	// Cutting top of EC
		if(!pFiducialCut(h_mom.Vect())) continue;
		h1_en_th_cut     -> Fill(    theta);
		h1_en_phi_cut    -> Fill(      phi);
		h2_en_th_phi_cut -> Fill(phi,theta);

		// ---------------------------------------
                // EC 10 cm cut
                h1_en_EC_x  -> Fill(ec_x[2]           );
                h1_en_EC_y  -> Fill(ec_y[2]           );
                h2_en_EC_xy -> Fill(ec_x[2], ec_y[2]  );

		// ---------------------------------------
		// Momentum cut
		h1_n_p     -> Fill(h_mom.Vect().Mag(),n_weight);
		if(h_mom.Vect().Mag()>2.34) continue;
		h1_n_p_cut -> Fill(h_mom.Vect().Mag(),n_weight);
		// ------------------------------------------------------------------------------

		Xb = Q2/(2.*Nu*Mn);

		u1 = h_mom.Vect().Unit();
		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;                
		p_over_q = h_mom.Rho()/q4.Rho();
		Mmiss    = (q4+n_pair-h_mom).Mag();

		if(Xb<1.1) continue;

		h1_Q2_Xb_cut_n -> Fill( Q2 , n_weight); 
		h2_thpq_npq    -> Fill(p_over_q,theta_pq);

		h2_thnq_pn    -> Fill(theta_pq,h_mom.Vect().Mag());
        	h2_thnq_thn   -> Fill(theta_pq,theta             );
        	h2_thnq_pmiss -> Fill(theta_pq,Pmiss             );
        	h2_thnq_Q2    -> Fill(theta_pq,Q2                );
        	h2_thnq_Nu    -> Fill(theta_pq,Nu                );
        	h2_thnq_Xb    -> Fill(theta_pq,Xb                );

		if((theta_pq < th_Nq_high )&&(theta_pq > th_Nq_low )&&(p_over_q > n_p_q_low)&&(p_over_q < n_p_q_high)){
			h1_Mmiss_n -> Fill(Mmiss,n_weight);
			h1_Pmiss_n -> Fill(Pmiss,n_weight);

			if((Pmiss>Pmiss_cut_lim[min_i])&&(Mmiss<Mmiss_cut_lim[min_j])&&(Pmiss<Pmiss_high)){
				h1_Emiss_n -> Fill(Emiss,n_weight);

				h1_pe_n -> Fill(e_mom.Rho(), n_weight );
				h1_te_n -> Fill(theta_e    , n_weight );
				h1_Q2_n -> Fill(Q2         , n_weight );
				h1_Nu_n -> Fill(Nu         , n_weight );
				h1_Xb_n -> Fill(Xb         , n_weight );

				h1_pm_n  -> Fill(Pmiss      , n_weight);
				h1_ph_n  -> Fill(h_mom.Rho(), n_weight);
				h1_th_n  -> Fill(theta      , n_weight);
				h1_thq_n -> Fill(theta_pq   , n_weight);

				if(p_points_to_EC_fid((q4+V4_Targ).Vect(),vtx_z[0],10.)){
					raw_nu_n += 1.;
					nu_een   += n_weight;
				}
			}
		}

	}
	cout << "**************************\n\n";

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Preparing plots
	TLegend ** leg_ec_cut = new TLegend*[4];
	for(int i = 0 ; i < 4 ; i++){
		leg_ec_cut[i] = new TLegend(0.13,0.7,0.3,0.85);
		leg_ec_cut[i] -> SetLineColor(0);
	}
	leg_ec_cut[0] -> AddEntry(h1_ep_EC_x     ,"p evts before cut");
	leg_ec_cut[0] -> AddEntry(h1_ep_EC_x_cut ,"p evts after  cut");
	leg_ec_cut[0] -> AddEntry(h1_en_EC_x     ,"n evts"           );
	leg_ec_cut[1] -> AddEntry(h1_ep_EC_y     ,"p evts before cut");
	leg_ec_cut[1] -> AddEntry(h1_ep_EC_y_cut ,"p evts after  cut");
	leg_ec_cut[1] -> AddEntry(h1_en_EC_y     ,"n evts"           );
	leg_ec_cut[2] -> AddEntry(h2_ep_EC_xy    ,"p evts before cut");
	leg_ec_cut[2] -> AddEntry(h2_ep_EC_xy_cut,"p evts after  cut");
	leg_ec_cut[3] -> AddEntry(h2_en_EC_xy    ,"n evts"           );
	leg_ec_cut[3] -> AddEntry(h2_ep_EC_xy_cut,"p evts after  cut");
	// ----
	h2_ep_EC_xy -> SetTitle("protons before and after 10 cm EC cut"  );
	h2_en_EC_xy -> SetTitle("protons and neutrons after 10 cm EC cut");
	// ----
	TLegend ** leg_th_phi_cut = new TLegend*[4];
	for(int i = 0 ; i < 4 ; i++){
		leg_th_phi_cut[i] = new TLegend(0.13,0.7,0.3,0.85);
		leg_th_phi_cut[i] -> SetLineColor(0);
	}
	leg_th_phi_cut[0] -> AddEntry(h1_en_th        ,"n evts before cut");
	leg_th_phi_cut[0] -> AddEntry(h1_en_th_cut    ,"n evts after  cut");
	leg_th_phi_cut[0] -> AddEntry(h1_ep_th        ,"p evts"           );
	leg_th_phi_cut[1] -> AddEntry(h1_en_phi       ,"n evts before cut");
	leg_th_phi_cut[1] -> AddEntry(h1_en_phi_cut   ,"n evts after  cut");
	leg_th_phi_cut[1] -> AddEntry(h1_ep_phi       ,"p evts"           );
	leg_th_phi_cut[2] -> AddEntry(h2_en_th_phi    ,"n evts before cut");
	leg_th_phi_cut[2] -> AddEntry(h2_en_th_phi_cut,"n evts after  cut");
	leg_th_phi_cut[3] -> AddEntry(h2_ep_th_phi    ,"p evts"           );
	leg_th_phi_cut[3] -> AddEntry(h2_en_th_phi_cut,"n evts after  cut");
	// ----
	TLegend * leg_p_cut = new TLegend(0.13,0.7,0.3,0.85);
	leg_p_cut -> SetLineColor(0);
	leg_p_cut -> AddEntry(h1_p_p    ,"p evts before cut");
	leg_p_cut -> AddEntry(h1_p_p_cut,"p evts after  cut");
	leg_p_cut -> AddEntry(h1_n_p    ,"n evts before cut");
	leg_p_cut -> AddEntry(h1_n_p_cut,"n evts after  cut");
	// ----
	TBox *b_thpq_ppq_cuts = new TBox(p_p_q_low,th_Nq_low,p_p_q_high,th_Nq_high);
	b_thpq_ppq_cuts -> SetFillStyle(0);
	b_thpq_ppq_cuts -> SetFillColor(0);
	b_thpq_ppq_cuts -> SetLineColor(2);
	b_thpq_ppq_cuts -> SetLineWidth(3);
	// ----
	TBox *b_thpq_npq_cuts = new TBox(n_p_q_low,th_Nq_low,n_p_q_high,th_Nq_high);
	b_thpq_npq_cuts -> SetFillStyle(0);
	b_thpq_npq_cuts -> SetFillColor(0);
	b_thpq_npq_cuts -> SetLineColor(2);
	b_thpq_npq_cuts -> SetLineWidth(3);
	// ----
	TLegend * leg_Q2 = new TLegend(0.6,0.7,0.85,0.85);
	leg_Q2 -> SetLineColor(0);
	leg_Q2 -> AddEntry(h1_Q2_Xb_cut_psm,"Smeared Protons");
	leg_Q2 -> AddEntry(h1_Q2_Xb_cut_n  ,"Neutrons"       );
	// ----
	TLegend * leg_m_p_miss = new TLegend(0.6,0.7,0.85,0.85);
	leg_m_p_miss -> SetLineColor(0);
	leg_m_p_miss -> AddEntry(h1_Mmiss_psm,"Smeared Protons");
	leg_m_p_miss -> AddEntry(h1_Mmiss_n  ,"Neutrons"       );
	// ----
	TLegend * leg_pm_fpn = new TLegend(0.6,0.7,0.85,0.85);
	leg_pm_fpn -> SetLineColor(0);
	leg_pm_fpn -> AddEntry(h1_Pm_psm   ,"Smeared p_{miss}");
	leg_pm_fpn -> AddEntry(h1_Pm_psm_fp,"False positive"  );
	leg_pm_fpn -> AddEntry(h1_Pm_psm_fn,"False negative"  );
	// ----
	double norm = (h1_Q2_Xb_cut_n->Integral())/(h1_Q2_Xb_cut_psm->Integral());
	h1_Q2_Xb_cut_psm -> Scale(norm);
	// ----
	norm = (h1_Mmiss_n->Integral())/(h1_Mmiss_psm->Integral());
	h1_Mmiss_psm -> Scale(norm);
	// ----
	norm = (h1_Pmiss_n->Integral())/(h1_Pmiss_psm->Integral());
	h1_Pmiss_psm -> Scale(norm);
	// ----
	norm = (h1_Emiss_n->Integral())/(h1_Emiss_psm->Integral());
	h1_Emiss_psm -> Scale(norm);
	// ----
	norm = (h1_pe_n->Integral())/(h1_pe_psm->Integral());
	h1_pe_psm -> Scale(norm);
	norm = (h1_te_n->Integral())/(h1_te_psm->Integral());
	h1_te_psm -> Scale(norm);
	norm = (h1_Q2_n->Integral())/(h1_Q2_psm->Integral());
	h1_Q2_psm -> Scale(norm);
	norm = (h1_Nu_n->Integral())/(h1_Nu_psm->Integral());
	h1_Nu_psm -> Scale(norm);
	norm = (h1_Xb_n->Integral())/(h1_Xb_psm->Integral());
	h1_Xb_psm -> Scale(norm);
	// ----
	norm = (h1_pm_n ->Integral())/(h1_pm_psm ->Integral());
	h1_pm_psm -> Scale(norm);
	norm = (h1_ph_n ->Integral())/(h1_ph_psm ->Integral());
	h1_ph_psm -> Scale(norm);
	norm = (h1_th_n ->Integral())/(h1_th_psm ->Integral());
	h1_th_psm -> Scale(norm);
	norm = (h1_thq_n->Integral())/(h1_thq_psm->Integral());
	h1_thq_psm -> Scale(norm);

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Editting plots for presentation
	pretty_TH1F(h1_Q2_Xb_cut_p   );
	pretty_TH1F(h1_Q2_Xb_cut_psm );
	pretty_TH1F(h1_Mmiss_p       );
	pretty_TH1F(h1_Mmiss_psm     );
	pretty_TH1F(h1_Pmiss_p       );
	pretty_TH1F(h1_Pmiss_psm     );
	pretty_TH1F(h1_Pm_psm        );
	pretty_TH1F(h1_Emiss_p       );
	pretty_TH1F(h1_Emiss_n       );
	pretty_TH1F(h1_pe_n          );
	pretty_TH1F(h1_te_n          );
	pretty_TH1F(h1_Q2_n          );
	pretty_TH1F(h1_Nu_n	     );
	pretty_TH1F(h1_Xb_n          );
	pretty_TH1F(h1_pm_psm        );
	pretty_TH1F(h1_ph_psm        );
	pretty_TH1F(h1_th_psm        );
	pretty_TH1F(h1_thq_psm       );

	pretty_TH2F(h2_thpq_ppq      );
	pretty_TH2F(h2_thpq_psmpq    );
	pretty_TH2F(h2_thpq_npq      );

	pretty_TH2F(h2_thnq_pn       );
        pretty_TH2F(h2_thnq_thn      );
        pretty_TH2F(h2_thnq_pmiss    );
        pretty_TH2F(h2_thnq_Q2       );
        pretty_TH2F(h2_thnq_Nu       );
        pretty_TH2F(h2_thnq_Xb       );

	for(int j = 0 ; j < no_mmiss_cuts ; j++){
		pretty_TGraphErrors(f_fal_pos[j]);
		pretty_TGraphErrors(f_fal_neg[j]);
		pretty_TGraphErrors(f_pos_neg[j]);
	}

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Plotting histograms
	gStyle->SetOptStat(0);

	TCanvas * c1 = new TCanvas("c1","",1300,800);
	c1 -> Divide(2,2);
	c1 -> cd(1);    h1_ep_EC_x  -> Draw();  h1_ep_EC_x_cut  -> Draw("same");        h1_en_EC_x -> Draw("same");     leg_ec_cut[0]-> Draw("same");
	c1 -> cd(2);    h1_ep_EC_y  -> Draw();  h1_ep_EC_y_cut  -> Draw("same");        h1_en_EC_y -> Draw("same");     leg_ec_cut[1]-> Draw("same");
	c1 -> cd(3);    h2_ep_EC_xy -> Draw();  h2_ep_EC_xy_cut -> Draw("same");                                        leg_ec_cut[2]-> Draw("same");
	c1 -> cd(4);    h2_en_EC_xy -> Draw();  h2_ep_EC_xy_cut -> Draw("same");                                        leg_ec_cut[3]-> Draw("same");
	c1->Modified(); c1->Update();

	TCanvas * c2 = new TCanvas("c2","",1300,800);
	c2 -> Divide(2,2);
	c2 -> cd(1);    h1_ep_th     -> Draw(); h1_en_th         -> Draw("same");       h1_en_th_cut  -> Draw("same");  leg_th_phi_cut[0]-> Draw("same");
	c2 -> cd(2);    h1_ep_phi    -> Draw(); h1_en_phi        -> Draw("same");       h1_en_phi_cut -> Draw("same");  leg_th_phi_cut[1]-> Draw("same");
	c2 -> cd(3);    h2_en_th_phi -> Draw(); h2_en_th_phi_cut -> Draw("same");                                       leg_th_phi_cut[2]-> Draw("same");
	c2 -> cd(4);    h2_ep_th_phi -> Draw(); h2_en_th_phi_cut -> Draw("same");                                       leg_th_phi_cut[3]-> Draw("same");
	c2->Modified(); c2->Update();

	TCanvas * c3 = new TCanvas("c3","",1300,800);
	c3 -> Divide(2,2);
	c3 -> cd(1);	h1_Q2_Xb_cut_p   -> Draw();
	c3 -> cd(2);    h1_Q2_Xb_cut_psm -> Draw();	h1_Q2_Xb_cut_n  -> Draw("same");	leg_Q2 -> Draw("same");
	c3->Modified(); c3->Update();

	TCanvas * c4 = new TCanvas("c4","",1300,800);
	c4 -> Divide(2,2);
	c4 -> cd(1);	h2_thpq_ppq   -> Draw("COLZ");	b_thpq_ppq_cuts -> Draw("same");
	c4 -> cd(3);    h2_thpq_psmpq -> Draw("COLZ");	b_thpq_npq_cuts -> Draw("same");
	c4 -> cd(4);    h2_thpq_npq   -> Draw("COLZ");	b_thpq_npq_cuts -> Draw("same");
	c4->Modified(); c4->Update();

	TCanvas * c5 = new TCanvas("c5","",1300,800);
	c5 -> Divide(2,2);
	c5 -> cd(1);	h1_Mmiss_p   -> Draw();
	c5 -> cd(2);    h1_Mmiss_psm -> Draw();	h1_Mmiss_n -> Draw("same");	leg_m_p_miss -> Draw("same");
	c5 -> cd(3);    h1_Pmiss_p   -> Draw();
	c5 -> cd(4);    h1_Pmiss_psm -> Draw(); h1_Pmiss_n -> Draw("same");	leg_m_p_miss -> Draw("same");
	c5->Modified(); c5->Update();

	TCanvas * c6 = new TCanvas("c6","",1300,800);
	c6 -> Divide(2,2);
	c6 -> cd(1);
	f_fal_pos[no_mmiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_mmiss_cuts-1 ; j++) f_fal_pos[j] -> Draw("samePE");
	//leg_false_pos -> Draw("same");
	c6 -> cd(2);
	f_fal_neg[no_mmiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_mmiss_cuts-1 ; j++) f_fal_neg[j] -> Draw("samePE");
	leg_false_neg -> Draw("same");
	c6 -> cd(3);
	f_pos_neg[no_mmiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_mmiss_cuts-1 ; j++) f_pos_neg[j] -> Draw("samePE");
	leg_pos_neg -> Draw("same");	l_f_pos_neg -> Draw("same");
	c6 -> cd(4);	h1_Pm_psm -> Draw();	h1_Pm_psm_fp -> Draw("same");	h1_Pm_psm_fn -> Draw("same");	leg_pm_fpn -> Draw("same");
	c6->Modified(); c6->Update();

	TCanvas * c7 = new TCanvas("c7","",1300,800);
	c7 -> Divide(2,2);
	c7 -> cd(1);	h1_Emiss_p -> Draw();
	c7 -> cd(2);    h1_Emiss_n -> Draw();	h1_Emiss_psm -> Draw("same");	leg_m_p_miss -> Draw("same");
	c7->Modified(); c7->Update();

	TCanvas * c8 = new TCanvas("c8","",1300,800);
	c8 -> Divide(3,2);
	c8 -> cd(1);	h1_pe_n -> Draw();	h1_pe_psm -> Draw("same");
	c8 -> cd(2);	h1_te_n -> Draw();	h1_te_psm -> Draw("same");
	c8 -> cd(3);	h1_Q2_n -> Draw();	h1_Q2_psm -> Draw("same");	leg_m_p_miss -> Draw("same");
	c8 -> cd(4);	h1_Nu_n -> Draw();	h1_Nu_psm -> Draw("same");
	c8 -> cd(5);	h1_Xb_n -> Draw();	h1_Xb_psm -> Draw("same");
	c8->Modified(); c8->Update();

	TCanvas * c9 = new TCanvas("c9","",1300,800);
	c9 -> Divide(2,2);
	c9 -> cd(1);	h1_pm_psm -> Draw();	h1_pm_n -> Draw("same");	leg_m_p_miss -> Draw("same");
	c9 -> cd(2);	h1_ph_psm -> Draw();	h1_ph_n -> Draw("same");
	c9 -> cd(3);	h1_th_psm -> Draw();	h1_th_n -> Draw("same");
	c9 -> cd(4);	h1_thq_psm-> Draw();	h1_thq_n-> Draw("same");
	c9->Modified(); c9->Update();

	TCanvas * c10 = new TCanvas("c10","",1300,800);
        c10 -> Divide(3,2);
       	c10 -> cd(1);	h2_thnq_pn    -> Draw("COLZ");
	c10 -> cd(2);	h2_thnq_thn   -> Draw("COLZ");
	c10 -> cd(3);	h2_thnq_pmiss -> Draw("COLZ");
	c10 -> cd(4);	h2_thnq_Q2    -> Draw("COLZ");
	c10 -> cd(5);	h2_thnq_Nu    -> Draw("COLZ");
	c10 -> cd(6);	h2_thnq_Xb    -> Draw("COLZ");
	c10->Modified(); c10->Update();

	// --------------------------
	TCanvas * c11 = new TCanvas("c11","",1300,800);
	c11 -> Divide(2,2);
	c11 -> cd(1);

	double g_model_ratio_x[2] = {N+Z-2,N+Z+2};
	double g_model_ratio_y[2];
        double g_model_ratio_e[2];
	prep_pink_band( nu_eep_een_GK05/nu_eep , nu_eep_een_AMT/nu_eep , nu_eep_een_MMD/nu_eep , g_model_ratio_y , g_model_ratio_e );
	TGraphErrors * g_model_ratio = new TGraphErrors(2,g_model_ratio_x,g_model_ratio_y,0,g_model_ratio_e);
	g_model_ratio -> SetFillStyle(3008);
	g_model_ratio -> SetFillColor(6);
	g_model_ratio -> GetYaxis() -> SetTitle("#frac{#(e,e'p)}{#(e,e'n)}");
	g_model_ratio -> GetXaxis() -> SetTitle("A");
	g_model_ratio -> SetTitle(tab_targ);
	g_model_ratio -> GetXaxis() -> SetRangeUser(g_model_ratio_x[0],g_model_ratio_x[1]);
	g_model_ratio -> SetMinimum(1.0);
	g_model_ratio -> SetMaximum(3.5);
	pretty_TGraphErrors(g_model_ratio);
	g_model_ratio -> Draw("AE3");

	double A_x = N+Z;
	double eep_een      = (nu_eep)/(nu_een);
	double e_sta_eep_een = sqrt( pow(sqrt(nu_eep)/nu_een,2) + pow(nu_eep*sqrt(nu_een)/pow(nu_een,2),2) );       // statistical error
	double e_def_eep_een = 0; //sqrt( pow(n_eff_p0*n_eff_m_e,2)+pow(n_eff_p0_e*n_eff_m,2) );                             // neutron detection eff error
	double e_tot_eep_een = sqrt(e_sta_eep_een*e_sta_eep_een+e_def_eep_een*e_def_eep_een);

	cout << "**************\n[#(e,e'p)]/[#(e,e'n)] = " << eep_een << " +/- " << e_tot_eep_een;
	cout << "\nstatistical error: " << e_sta_eep_een;
	cout << "\nn det eff   error: " << e_def_eep_een << endl;
	cout << "**************" << endl;

	cout << "Number of protons  contributing to this ratio: " << raw_nu_p << endl;
	cout << "Number of neutrons contributing to this ratio: " << raw_nu_n << endl << endl;

	TGraphErrors * g_data_ratio = new TGraphErrors(1,&A_x,&eep_een,0,&e_tot_eep_een);
	g_data_ratio -> SetMarkerStyle(20);
	g_data_ratio -> Draw("sameP");

	c11->Modified(); c11->Update();

	// --------------------------
	// Printing plots to pdf
	c3  ->Print("./Results/plots_SRC_"+tab_targ+".pdf(");
	c4  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c5  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c6  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c7  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c8  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c9  ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c10 ->Print("./Results/plots_SRC_"+tab_targ+".pdf" );
	c11 ->Print("./Results/plots_SRC_"+tab_targ+".pdf)");

	// --------------------------
	// Creating output table
	ofstream output_table;
	output_table.open("./Results/table_SRC_"+tab_targ+".dat");
	output_table << tab_targ      << endl;
	output_table << eep_een       << endl;
	output_table << e_tot_eep_een << endl;
	output_table.close();

	myapp->Run();
	return 0;
}
// ====================================================================================================================================================
void load_cross_section_maps(){
        // Opening input file
        TFile * TF2_sigma = new TFile("./parameters/TF2_cross_section_maps.root");
        sigma_p_GK05 = (TF2*) TF2_sigma->Get("sigma_p_GK05")->Clone();
        sigma_n_GK05 = (TF2*) TF2_sigma->Get("sigma_n_GK05")->Clone();
        sigma_p_AMT  = (TF2*) TF2_sigma->Get("sigma_p_AMT" )->Clone();
        sigma_n_AMT  = (TF2*) TF2_sigma->Get("sigma_n_AMT" )->Clone();
        sigma_p_MMD  = (TF2*) TF2_sigma->Get("sigma_p_MMD" )->Clone();
        sigma_n_MMD  = (TF2*) TF2_sigma->Get("sigma_n_MMD" )->Clone();
        TF2_sigma -> Close();
	cout << "Loaded TF2 functions for single-nucleon cross-section calculations" << endl;
}
// ====================================================================================================================================================
void prep_pink_band( double GK05_res , double AMT_res , double MMD_res , double * y_ctr_array , double * y_err_array ){
        double max_temp = TMath::Max( GK05_res , AMT_res  );
        double max_perm = TMath::Max( MMD_res  , max_temp );
        double min_temp = TMath::Min( GK05_res , AMT_res  );
        double min_perm = TMath::Min( MMD_res  , min_temp );

        double y_ctr = (max_perm + min_perm )/2.;
        double y_err = (max_perm - min_perm )/2.;

        for(int i = 0 ; i < 2 ; i++){
                y_ctr_array[i] = y_ctr;
                y_err_array[i] = y_err;
        }
}
// ====================================================================================================================================================
void pretty_TH2F(TH2F * h2){
	h2 -> GetXaxis() -> SetLabelSize(0.08);
	h2 -> GetXaxis() -> SetNdivisions(505);
	h2 -> GetXaxis() -> SetTitleSize(0.08);
	h2 -> GetXaxis() -> CenterTitle();

	h2 -> GetYaxis() -> SetLabelSize(0.08);
	h2 -> GetYaxis() -> SetNdivisions(505);
	h2 -> GetYaxis() -> SetTitleSize(0.08);
	h2 -> GetYaxis() -> CenterTitle();

	gStyle->SetPadBottomMargin(0.18);
	gStyle->SetPadLeftMargin  (0.18);
}
// ====================================================================================================================================================
void pretty_TH1F(TH1F * h1){
	h1 -> GetXaxis() -> SetLabelSize(0.08);
	h1 -> GetXaxis() -> SetNdivisions(505);
	h1 -> GetXaxis() -> SetTitleSize(0.08);
	h1 -> GetXaxis() -> CenterTitle();

	h1 -> GetYaxis() -> SetLabelSize(0.08);
	h1 -> GetYaxis() -> SetNdivisions(505);
	h1 -> GetYaxis() -> SetTitleSize(0.08);
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetTitleOffset(1.2);

	gStyle->SetPadBottomMargin(0.18);
	gStyle->SetPadLeftMargin  (0.18);
}
// ====================================================================================================================================================
void pretty_TGraphErrors(TGraphErrors * gP){
	gP -> GetXaxis() -> SetLabelSize(0.08);
	gP -> GetXaxis() -> SetNdivisions(505);
	gP -> GetXaxis() -> SetTitleSize(0.08);
	gP -> GetXaxis() -> CenterTitle();

	gP -> GetYaxis() -> SetLabelSize(0.08);
	gP -> GetYaxis() -> SetNdivisions(505);
	gP -> GetYaxis() -> SetTitleSize(0.08);
	gP -> GetYaxis() -> CenterTitle();

	gStyle->SetPadBottomMargin(0.18);
	gStyle->SetPadLeftMargin  (0.18);
}
// ====================================================================================================================================================
float Coulomb_Corr_MeV(int Z, int A){
	// See Effective Momentum Approximation (EMA) formalism
	// https://www.jlab.org/indico/event/206/session/27/contribution/12/material/slides/0.pdf

	const float alpha  = 0.0072973525664;
	const float hbar_c = 197.327; // MeV*fm

	float R0_fm = 1.1*pow((float)A,1/3.) + 0.86*pow((float)A,-1/3.); // nuclear radius [fm]
	float DV = 3*((float)Z-1)*alpha*hbar_c/(2.*R0_fm);
	float DE = 0.775*DV;

	return DE;
}
// ====================================================================================================================================================
TLorentzVector P4_Coulomb_Corrected(TLorentzVector P4_Uncorrected, int Z, int A){
	float Correction_GeV = Coulomb_Corr_MeV( Z , A )/1000.;
	float E_corr  = P4_Uncorrected.E() - Correction_GeV;
	float New_Mag = sqrt(E_corr*E_corr-Mp*Mp);

	TVector3 unit = P4_Uncorrected.Vect().Unit();
	unit.SetMag(New_Mag);

	P4_Uncorrected.SetXYZM(unit.X(),unit.Y(),unit.Z(),Mp);
	return P4_Uncorrected;
}
// ====================================================================================================================================================
void read_n_mom_res_params(){
	ifstream n_mom_res_params;
	n_mom_res_params.open("./parameters/EC_p_res.dat");
	n_mom_res_params >> ec_time_res;
	n_mom_res_params >> ec_extr_res;
	n_mom_res_params >> ec_extr_res_err;
	n_mom_res_params.close();
	cout<<"Loaded parameters for neutron momentum resolution ("<<ec_time_res<<", "<<ec_extr_res<<" +/- "<<ec_extr_res_err<<")"<<endl;
}
// ====================================================================================================================================================
void read_n_det_eff_params(){
        TFile * neutron_detection_efficiency = new TFile("./parameters/neut_det_eff.root");
        n_det_eff     = (TGraph*) neutron_detection_efficiency -> Get("neut_det_eff"   )  -> Clone();
        n_det_eff_err = (TGraph*) neutron_detection_efficiency -> Get("neut_det_eff_err") -> Clone();
        neutron_detection_efficiency -> Close();
        cout<<"Loaded parameters for neut. det. efficiency" << endl;
}
// ====================================================================================================================================================
double n_efficiency     ( double mom_GeV ){
        if      ((mom_GeV>=0.5)&&(mom_GeV<=2.0)) return n_det_eff -> Eval(mom_GeV);
        else if                 ( mom_GeV> 2.0 ) return n_det_eff -> Eval( 2.0   );
        else{cout << "Neutron detection efficiency not determined in the given momentum range" << endl; exit(-1);}
}
// ====================================================================================================================================================
double n_efficiency_err ( double mom_GeV ){
        if      ((mom_GeV>=0.5)&&(mom_GeV<=2.0)) return n_det_eff_err -> Eval(mom_GeV);
        else if                 ( mom_GeV> 2.0 ) return n_det_eff_err -> Eval( 2.0   );
        else{cout << "Neutron detection efficiency not determined in the given momentum range" << endl; exit(-1);}
}
// ====================================================================================================================================================
double f_delta_p(double p){
	// Returns value of delta_p as a function of p.
	// This function is used for neutron momentum resolution.
	double dt = ec_time_res;
	const double d      = 516.9152;       //cm, distance from interaction vtx to EC
	const double Mn     = 0.939565378;    // Neutron mass in GeV
	const double c_cm_s = 2.99792458E+10; // cm/s, speed of light
	const double ns_2_s = 1E-09;

	double val = p*p/(Mn*Mn)*sqrt(Mn*Mn+p*p)*dt/d*c_cm_s*ns_2_s;
	return sqrt(val*val/p/p+ec_extr_res*ec_extr_res)*p;
}
// ====================================================================================================================================================
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;								// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;								// Missing energy
}
// ====================================================================================================================================================
bool p_points_to_EC_fid(TVector3 pm, double vtx_z, double dist){
	double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;

	const double ec_theta =  65.*pi/180.; // degrees
	const double z0       =563.1; // cm

	// Using law of sines to determine distance from vertex to EC
	pm_theta    = pm.Theta();
	third_theta = (pi - ec_theta - pm_theta);

	r = (z0-vtx_z)*sin(ec_theta)/sin(third_theta);

	double Xmiss = r*sin(pm.Theta())*cos(pm.Phi());
	double Ymiss = r*sin(pm.Theta())*sin(pm.Phi());

	pm_phi = 180./pi*pm.Phi();
	if (pm_phi < -30.) pm_phi += 360.;
	int sec = (int)(pm_phi+30)/60;
	if (sec>5) sec = 5;
	if (sec<0) sec = 0;

	angle = pi/180.*60.*(sec);

	x_rot = Xmiss*cos(angle) + Ymiss*sin(angle);
	y_rot = Xmiss*sin(angle) - Ymiss*cos(angle);

	d_intercept = dist/cos(atan(1.73));

	return((x_rot<390.-dist)&&(x_rot>1.73*abs(y_rot)+55.+d_intercept));

}
// ====================================================================================================================================================
bool read_p_fid_params(){
	//Parameters for 4 GeV proton's Fiducial Cut Rustam Niyazov
	//"http://www.physics.odu.edu/~rust/clas/fidp.html"

	ifstream param_file;
	param_file.open("./parameters/PFID_4461_2250.dat");

	for(int i = 0 ; i < 6 ; i++){
		for(int j = 0 ; j < 6 ; j++){
			param_file >> fgPar_Pfidft1l[i][j];
			param_file >> fgPar_Pfidft1r[i][j];
			param_file >> fgPar_Pfidft2l[i][j];
			param_file >> fgPar_Pfidft2r[i][j];
			param_file >> fgPar_Pfidbt1l[i][j];
			param_file >> fgPar_Pfidbt1r[i][j];
			param_file >> fgPar_Pfidbt2l[i][j];
			param_file >> fgPar_Pfidbt2r[i][j];
			param_file >> fgPar_Pfidbl  [i][j];
			param_file >> fgPar_Pfidbr  [i][j];
		}
	}      

	param_file.close();
	return true;
}
// ====================================================================================================================================================
bool pFiducialCut(TVector3 momentum){
	//Positive Hadron Fiducial Cut
	//Check out "http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html"
	Bool_t status = kTRUE;

	Float_t theta = momentum.Theta()*180/pi;
	Float_t phi   = momentum.Phi()  *180/pi;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();

	// ----------------------------------------------------------------------------------------------------------------
	Float_t parfidl [3];    for(Int_t i=0; i<3; i++){parfidl [i]=0;}
	Float_t parfidr [3];    for(Int_t i=0; i<3; i++){parfidr [i]=0;}
	Float_t parfidbl[2];    for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
	Float_t parfidbr[2];    for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
	Float_t cphil =0;       Float_t cphir =0;
	Float_t phi45l=0;       Float_t phi45r=0;
	Float_t phi60l=0;       Float_t phi60r=0;

	Float_t theta_min = 11;
	Float_t theta_max =140;

	bool Forward=kFALSE;            // defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
	Int_t   thetab    =45;          // this variable defines the edge point for Forward<->Backward regions
	Float_t p1        =0.575;       // last bin momentum for region p<0.6 GeV/c
	if(p<0.2) p=0.2;                // momentum less than 0.2 GeV/c, use 0.2 GeV/c
	if(p>4.4) p=4.4;                // momentum greater than 4.4 GeV/c, use 4.4 GeV/c
	// ----------------------------------------------------------
	//Get parametrized values of theta_max for p<0.6 GeV/c region
	if(p<0.6){theta_max=fgPar_Pfidbl[sector][4]+fgPar_Pfidbl[sector][5]*p;}
	//Get parametrized values of theta_max for p>0.6 GeV/c region
	else{theta_max=fgPar_Pfidbl[sector][4]+fgPar_Pfidbl[sector][5]*p1;}

	//Get the momentum dependent parameters for Forward Region (theta <45 deg)   
	Forward=kTRUE;
	if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
		//parameters for hyperbolic function
		for (Int_t i=0; i<3; i++){
			Int_t j=2*i;
			parfidl[i]=fgPar_Pfidft1l[sector][j]+fgPar_Pfidft1l[sector][j+1]/p;
			parfidr[i]=fgPar_Pfidft1r[sector][j]+fgPar_Pfidft1r[sector][j+1]/p;
		}
	}
	else{//forward2 defines  regions of momenta and p>0.6 GeV/c
		for (Int_t i=0; i<3; i++){
			Int_t j=2*i;
			parfidl[i]=fgPar_Pfidft2l[sector][j]+fgPar_Pfidft2l[sector][j+1]/p;
			parfidr[i]=fgPar_Pfidft2r[sector][j]+fgPar_Pfidft2r[sector][j+1]/p;
		}
	}
	// ----------------------------------------------------------
	phi45l=parfidl [0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
	phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
	if(theta>thetab){//backward region defined by theta >45 deg. 
		if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
		if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

		//Get the momentum dependent parameters for Backward Region

		Forward=kFALSE;
		if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
			//parameters for quadratic function
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_Pfidbt1l[sector][j]+fgPar_Pfidbt1l[sector][j+1]/p;
				parfidr[i]=fgPar_Pfidbt1r[sector][j]+fgPar_Pfidbt1r[sector][j+1]/p;
			}
			//these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
			for (Int_t i=0; i<2; i++){
				Int_t j=2*i;
				parfidbl[i]=fgPar_Pfidbl[sector][j]+fgPar_Pfidbl[sector][j+1]/p;
				parfidbr[i]=fgPar_Pfidbr[sector][j]+fgPar_Pfidbr[sector][j+1]/p;
			}
		}
		else{//backward2 defines  regions of momenta p>0.6 GeV/c
			//parameters for quadratic function
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_Pfidbt2l[sector][j]+fgPar_Pfidbt2l[sector][j+1]/p;
				parfidr[i]=fgPar_Pfidbt2r[sector][j]+fgPar_Pfidbt2r[sector][j+1]/p;
			}
			//these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
			for (Int_t i=0; i<2; i++){
				Int_t j=2*i;
				parfidbl[i]=fgPar_Pfidbl[sector][j]+fgPar_Pfidbl[sector][j+1]/p1;
				parfidbr[i]=fgPar_Pfidbr[sector][j]+fgPar_Pfidbr[sector][j+1]/p1;
			}
		}
	}
	// -------------------------------------------------------------------------------------------------------------------
	if(Forward){//Forward region
		if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.   
		cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
		cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
	}
	// -------------------------------------------------------------------------------------------------------------------
	else{//Backward region
		phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
		phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

		if(theta<60){
			cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
			cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
		}
		Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0; 
		//dr and er are theta_flat and phi_edge parameters for phi>0;  
		dl=parfidbl[0];el=parfidbl[1];
		dr=parfidbr[0];er=parfidbr[1];

		if(theta>45&&theta<60){ //BackwardA region
			//try to match parametrized values from Forward region to Backward region parameters
			if(cphil>phi45l)cphil=phi45l;
			if(cphir<phi45r)cphir=phi45r;
		}
		//BackwardB region & phi<0
		else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant 
		else if(theta>dl&&theta<=theta_max){
			cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line 
		else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
		//BackwardB region & phi>0
		if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant 
		else if(theta>dr&&theta<=theta_max){
			cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line 
		else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
	}//Backward Region
	// -------------------------------------------------------------------------------------------------------------------
	if(phi<0) status=(phi>cphil); //check the constrains 
	else if(phi>=0) {status=(phi<cphir);
	}

	if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min
	if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

	//p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for 
	//some range of momentum, where edge does not look good.
	bool s1s4 =(theta<11.7&&(sector==0||sector==3));
	bool s5   =(theta<12.2&& sector==4);
	bool s6   =(theta<11.4&& sector==5);
	if( p>=0.6 && p<1.5 && (s1s4||s5||s6) ) status=kFALSE;
	// ----------------------------------------------------------------------------------------------------------------
	return status;
}

