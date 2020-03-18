// ==================================
// Code for the 4He/3He/12C analysis
// QE Mean-Field nucleons
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

const double BE_2H   = 0.00222;			// 2H  binding energy in GeV
const double BE_3H   = 0.00848;			// 3H  binding energy in GeV
const double BE_3He  = 0.00772;			// 3He binding energy in GeV
const double BE_4He  = 0.02830;			// 4He binding energy in GeV
const double BE_11B  = 0.07620;			// 11B binding energy in GeV
const double BE_11C  = 0.07344;			// 11C binding energy in GeV
const double BE_12C  = 0.09215;			// 12C binding energy in GeV

const double M2H     =   Mp +   Mn -BE_2H ;
const double M3H     =   Mp + 2*Mn -BE_3H ;
const double M3He    = 2*Mp +   Mn -BE_3He;	// 3He mass in GeV
const double M4He    = 2*Mp + 2*Mn -BE_4He;	// 4He mass in GeV
const double M11B    = 5*Mp + 6*Mn -BE_11B;
const double M11C    = 6*Mp + 5*Mn -BE_11C;
const double M12C    = 6*Mp + 6*Mn -BE_12C;	// 12C mass in GeV
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
void pretty_TH2F(TH2F * h2);								// Makes TH2F's look nice
void pretty_TH1F(TH1F * h1);								// Makes TH1F's look nice
void pretty_TGraphErrors(TGraphErrors * gP);						// Makes TGraphError's look nice
float Coulomb_Corr_MeV(int Z, int A);                                                   // Calculates proton Coulomb Energy Correction
TLorentzVector P4_Coulomb_Corrected(TLorentzVector P4_Uncorrected, int Z, int A);       // Applies proton Coulomb Energy Correction to proton momentum
double y_scale_func(TLorentzVector y_q, double y_Mtar, double y_Mres, double y_Mnuc);	// Returns value of y variable in units of GeV
void read_n_mom_res_params();								// Reads in parameters needed for neutron momentum resolution
void read_n_det_eff_params();								// Reads in parameters needed for neutron detection effocoency
double n_efficiency     ( double mom_GeV );						// Returns neutron detection efficiency
double n_efficiency_err ( double mom_GeV );						// Returns neutron detection efficiency error
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
	const double z0       =563.1; // cm
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
	double  Pmiss, Emiss   , p_ecx, p_ecy, p_smeared, pp_min_psm, Pmiss_sm, Emiss_sm, Mres, 
		y_sc , theta_pq, Pmiss_par, Xb, r_ec, pm_theta, third_theta, theta, phi, phi_e, theta_e,
		theta_e_calc, theta_h_calc, p_over_q, Z, N, n_weight;
	// --------------------
	// Variables for cross section ratios
	double nu_eep = 0;
	double nu_een = 0;
	double raw_nu_p = 0;
	double raw_nu_n = 0;
	const int bins_eepn = 3;
	double p_bins_eepn[bins_eepn+1];
	double p_ctr_eepn [bins_eepn  ];
	double p_err_eepn [bins_eepn  ];
	double nu_eep_p   [bins_eepn  ] = {0};
	double nu_een_p   [bins_eepn  ] = {0};
	for(int i = 0 ; i < bins_eepn+1 ; i++) {p_bins_eepn[i] = 1.6 + (float)i*(2.1 - 1.6)/(float)bins_eepn;}
	for(int i = 0 ; i < bins_eepn   ; i++) {p_ctr_eepn [i] = (p_bins_eepn[i]+p_bins_eepn[i+1])/2.;	p_err_eepn[i] = (p_bins_eepn[i+1]-p_bins_eepn[i])/2.;}
	// --------------------
	// Variables for theoretical cross section ratios
	double nu_eep_een_GK05 = 0;
	double nu_eep_een_AMT  = 0;
	double nu_eep_een_MMD  = 0;
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

	// QE cuts
	double Pmiss_cut = 0.25; // GeV
	double Emiss_cut = 0.08; // GeV 
	// --------------------------------------------------
	// electron kinematical cuts
	double y_scl_min = -0.08; // GeV
	double y_scl_max =  0.20; // GeV
	double Nu_min    =  0.90; // GeV
	double Nu_max    =  1.44; // GeV
	double Q2_min    =  1.20; // GeV^2
	double Q2_max    =  2.50; // GeV^2
	double Th_pq_max =  7.00; // deg	
	// --------------------------------------------------
	// Optimization of QE cuts
	const int no_emiss_cuts = 7;
	float Emiss_cut_lim[no_emiss_cuts] = {0};
	for(int i = 0 ; i < no_emiss_cuts ; i++) Emiss_cut_lim[i] = 0.04+0.05*(float)i;

	const int no_pmiss_cuts = 14;
	float Pmiss_cut_lim[no_pmiss_cuts] = {0};
	for(int i = 0 ; i < no_pmiss_cuts ; i++) Pmiss_cut_lim[i] = 0.09+0.021*(float)i;

	float fal_denom  [no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_pos_num[no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_pos    [no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_pos_e  [no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_neg_num[no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_neg    [no_emiss_cuts][no_pmiss_cuts] = {{0}};
	float fal_neg_e  [no_emiss_cuts][no_pmiss_cuts] = {{0}};

	// ----------------------------------------------------------
	//Define Histograms
	TH2F * h2_p_Pm_Em   = new TH2F("h2_p_Pm_Em"  ,"(a) protons (acc cuts);P_{miss} [GeV];E_{miss} [GeV]"        , 80 , 0.0 , 0.5 , 80 , 0.0 , 0.5 );
	TH2F * h2_n_Pm_Em   = new TH2F("h2_n_Pm_Em"  ,"(b) neutrons (acc cuts);P_{miss} [GeV];E_{miss} [GeV]"       , 80 , 0.0 , 0.5 , 80 , 0.0 , 0.5 );
	TH2F * h2_psm_Pm_Em = new TH2F("h2_psm_Pm_Em","(c) smeared protons (acc cuts);P_{miss} [GeV];E_{miss} [GeV]", 80 , 0.0 , 0.5 , 80 , 0.0 , 0.5 );

	TH2F * h2_p_Pm_Em_cut1   = new TH2F("h2_p_Pm_Em_cut1"  ,"(c) protons (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]"        ,80,0,.5,80,0,.5);
	TH2F * h2_n_Pm_Em_cut1   = new TH2F("h2_n_Pm_Em_cut1"  ,"(b) neutrons (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]"       ,80,0,.5,80,0,.5);
	TH2F * h2_psm_Pm_Em_cut1 = new TH2F("h2_psm_Pm_Em_cut1","(a) smeared protons (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]",80,0,.5,80,0,.5);
	TH2F * h2_psm_Pm_Em_cut2 = new TH2F("h2_psm_Pm_Em_cut2","(d) smeared protons (acc + p/E_{miss} + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]"
			, 80 , 0.0 , 0.5 , 80 , 0.0 , 0.5 );

	TH2F * h2_p_Nu_y          = new TH2F("h2_p_Nu_y"          ,"protons (acc cuts);#omega [GeV];y [GeV];Counts"             ,200,0.3,3.5,200,-0.5,0.7);
	TH2F * h2_p_Nu_y_QEcut    = new TH2F("h2_p_Nu_y_QEcut"    ,"protons (acc + p/E_{miss} cuts);#omega [GeV];y [GeV];Counts",200,0.3,3.5,200,-0.5,0.7);
	TH2F * h2_p_Nu_y_QEcut2   = new TH2F("h2_p_Nu_y_QEcut2"   ,"protons (acc + outside p/E_{miss} cuts);#omega [GeV];y [GeV];Counts",200,0.3,3.5,200,-0.5,0.7);

	TH2F * h2_p_Q2_thpq       = new TH2F("h2_p_Q2_thpq"       ,"protons (acc cuts);Q^{2} [GeV^{2}];#theta_{pq} [deg];Counts"              ,200,0,5.0,200,0,60);
	TH2F * h2_p_Q2_thpq_QEcut = new TH2F("h2_p_Q2_thpq_QEcut" ,"protons (acc + p/E_{miss} cuts);Q^{2} [GeV^{2}];#theta_{pq} [deg];Counts" ,200,0,5.0,200,0,60);
	TH2F * h2_p_Q2_thpq_QEcut2= new TH2F("h2_p_Q2_thpq_QEcut2","protons (acc + outside p/E_{miss} cuts);Q^{2} [GeV^{2}];#theta_{pq} [deg];Counts" ,200,0,5.0,200,0,60);

	TH1F * h1_y_cut2          = new TH1F("h1_y_cut2"     ,"y dist (acc + p/E_{miss} + e kinem cuts);y [GeV];Counts"                    ,100,-0.2, 0.3);
	TH1F * h1_Nu_cut2         = new TH1F("h1_Nu_cut2"    ,"#omega dist (acc + p/E_{miss} + e kinem cuts);#omega [GeV];Counts"          ,100, 0.8, 1.6);
	TH1F * h1_th_pq_cut2      = new TH1F("h1_th_pq_cut2" ,"#theta_{pq} dist (acc + p/E_{miss} + e kinem cuts);#theta_{pq} [deg];Counts",100, 0.0, 8.0);
	TH1F * h1_Q2_cut2         = new TH1F("h1_Q2_cut2"    ,"Q^{2} dist (acc + p/E_{miss} + e kinem cuts);Q^{2} [GeV^{2}];Counts"        ,100, 0.0, 4.0);
	// ----------------------------------
	// Acceptance matching: 10 cm EC cut
	TH2F * h2_en_EC_xy       = new TH2F("h2_en_EC_xy"      ,";EC_{x} [cm];EC_{y} [cm];Counts"      , 200 , -500, 500 , 200 , -500, 500 );
	TH2F * h2_ep_EC_xy       = new TH2F("h2_ep_EC_xy"      ,";EC_{x} [cm];EC_{y} [cm];Counts"      , 200 , -500, 500 , 200 , -500, 500 );
	TH2F * h2_ep_EC_xy_cut   = new TH2F("h2_ep_EC_xy_cut"  ,";EC_{x} [cm];EC_{y} [cm];Counts"      , 200 , -500, 500 , 200 , -500, 500 );

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
	TH2F * h2_ep_th_phi      = new TH2F("h2_ep_th_phi"      ,";#phi [deg];#theta [deg];Counts"    ,100,-100.,390.,100,0.,70.);
	TH2F * h2_en_th_phi      = new TH2F("h2_en_th_phi"      ,";#phi [deg];#theta [deg];Counts"    ,100,-100.,390.,100,0.,70.);
	TH2F * h2_en_th_phi_cut  = new TH2F("h2_en_th_phi_cut"  ,";#phi [deg];#theta [deg];Counts"    ,100,-100.,390.,100,0.,70.);

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
	// Smearing protons to simulate neutrons
	TH1F * h1_pp_mi_psm  = new TH1F("h1_pp_mi_psm"      ,";p_{p} - p_{smeared} [GeV];Counts" ,100,  -1.,1.  );

	TH2F * h2_p_pm_em    = new TH2F("h2_p_pm_em"    ,"protons zoomed out (acc cuts);P_{miss} [GeV];E_{miss} [GeV]"                  ,100,0.,0.5,100,-.6,.7);
	TH2F * h2_n_pm_em    = new TH2F("h2_n_pm_em"    ,"neutrons zoomed out (acc cuts);P_{miss} [GeV];E_{miss} [GeV]"                 ,100,0.,0.5,100,-.6,.7);
	TH2F * h2_psm_pm_em  = new TH2F("h2_psm_pm_em"  ,"smeared protons zoomed out (acc cuts);P_{miss} [GeV];E_{miss} [GeV]"          ,100,0.,0.5,100,-.6,.7);

	TH2F * h2_p_pm_em_1  = new TH2F("h2_p_pm_em_1"  ,"protons zoomed out (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]"        ,100,0.,0.5,100,-.6,.7);
	TH2F * h2_n_pm_em_1  = new TH2F("h2_n_pm_em_1"  ,"neutrons zoomed out (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]"       ,100,0.,0.5,100,-.6,.7);
	TH2F * h2_psm_pm_em_1= new TH2F("h2_psm_pm_em_1","smeared protons zoomed out (acc + e kinem cuts);P_{miss} [GeV];E_{miss} [GeV]",100,0.,0.5,100,-.6,.7);

	// Smeared - unsmeared correlations
	TH2F * h2_pm_pmsm = new TH2F("h2_pm_pmsm","(acc + e kinem cuts) smeared vs. unsmeared;P_{miss}^{unsmeared} [GeV];P_{miss}^{smeared} [GeV]",100, 0 ,1.8,100, 0,1.8);
	TH2F * h2_em_emsm = new TH2F("h2_em_emsm","(acc + e kinem cuts) smeared vs. unsmeared;E_{miss}^{unsmeared} [GeV];E_{miss}^{smeared} [GeV]",100,-.4,1.1,100,-1,1.2);
	TH2F * h2_ph_phsm = new TH2F("h2_ph_phsm","(acc + e kinem cuts) smeared vs. unsmeared;P_{p}^{unsmeared} [GeV];P_{p}^{smeared} [GeV]"      ,100, 0 ,2.4,100, 0,2.4);
	TH2F * h2_ph_pmsm = new TH2F("h2_ph_pmsm","(acc + e kinem cuts) smeared vs. unsmeared;P_{p}^{unsmeared} [GeV];P_{miss}^{smeared} [GeV]"   ,100, 0 ,2.4,100, 0,1.8);

	TH2F * h2_pm_pmsm_pm_diff = new TH2F("h2_pm_pmsm_pm_diff",
			"(acc + e kinem cuts) unsmeared vs. smeared-unsmeared;P_{miss}^{smeared}-P_{miss}^{unsmeared} [GeV];P_{miss}^{unsmeared} [GeV]",
			100,-0.5,1.2,100,  0,1.8);
	TH2F * h2_em_emsm_em_diff = new TH2F("h2_em_emsm_em_diff",
			"(acc + e kinem cuts) unsmeared vs. smeared-unsmeared;E_{miss}^{smeared}-E_{miss}^{unsmeared} [GeV];E_{miss}^{unsmeared} [GeV]",
			100,-1.0,1.0,100,-.4,1.1);
	TH2F * h2_ph_phsm_ph_diff = new TH2F("h2_ph_phsm_ph_diff",
			"(acc + e kinem cuts) unsmeared vs. smeared-unsmeared;P_{p}^{smeared}-P_{p}^{unsmeared} [GeV];P_{p}^{unsmeared} [GeV]"         ,
			100,-1.4,0.8,100,  0,2.4);
	TH2F * h2_ph_pmsm_ph_diff = new TH2F("h2_ph_pmsm_ph_diff",
			"(acc + e kinem cuts) unsmeared vs. smeared-unsmeared;P_{miss}^{smeared}-P_{p}^{unsmeared} [GeV];P_{p}^{unsmeared} [GeV]"      ,
			100,-2.5,1.5,100,  0,2.4);

	TH2F * h2_but_pm_em_0 = new TH2F("h2_but_pm_em_0",
			"(acc cuts);E_{miss}^{smeared}-E_{miss}^{unsmeared} [GeV];P_{miss}^{smeared}-P_{miss}^{unsmeared} [GeV]"          ,200,-0.8,1.5,200,-1,1.5);
	TH2F * h2_but_pm_em_1 = new TH2F("h2_but_pm_em_1",
			"(acc + e kinem cuts);E_{miss}^{smeared}-E_{miss}^{unsmeared} [GeV];P_{miss}^{smeared}-P_{miss}^{unsmeared} [GeV]",200,-0.8,1.5,200,-1,1.5);
	TH2F * h2_but_pm_ph   = new TH2F("h2_but_pm_ph",
			"(acc + e kinem cuts);P_{p}^{smeared}-P_{p}^{unsmeared} [GeV];P_{miss}^{smeared}-P_{miss}^{unsmeared} [GeV]",200,-1.5,1,200,-1,1.5);
	TH2F * h2_but_em_ph   = new TH2F("h2_but_em_ph",
			"(acc + e kinem cuts);P_{p}^{smeared}-P_{p}^{unsmeared} [GeV];E_{miss}^{smeared}-E_{miss}^{unsmeared} [GeV]",200,-1.5,1,200,-0.8,1.5);

	// ----------------------------------
	// Comparing smeared and un-smeared protons
	TH2F * h2_thp_y  = new TH2F("h2_thp_y"  ,"protons (p/E_{miss} cuts);#theta_{p} [deg];y [GeV];Counts" , 200 , 10, 80 , 200 , -0.5 , 0.6 );
	TH1F * h1_y_ec   = new TH1F("h1_y_ec"   ,"protons (p/E_{miss} cuts);y [GeV];Counts"                                 , 200 , -0.4 , 0.4 );
	TH1F * h1_y_full = new TH1F("h1_y_full" ,"protons (p/E_{miss} cuts);y [GeV];Counts"                                 , 200 , -0.4 , 0.4 );
	h1_y_ec -> SetLineColor(2);

	TH2F * h2_thp_pmpar  = new TH2F("h2_thp_pmpar"  ,"protons (p/E_{miss} cuts);#theta_{p} [deg];p^{||}_{miss} [GeV];Counts" , 200 , 10, 80 , 200 , -0.3 , 0.3 );
	TH1F * h1_pmpar_ec   = new TH1F("h1_pmpar_ec"   ,"protons (p/E_{miss} cuts);p^{||}_{miss} [GeV];Counts"                                 , 200 , -0.3 , 0.3 );
	TH1F * h1_pmpar_full = new TH1F("h1_pmpar_full" ,"protons (p/E_{miss} cuts);p^{||}_{miss} [GeV];Counts"                                 , 200 , -0.3 , 0.3 );
	h1_pmpar_ec -> SetLineColor(2);

	TH2F * h2_y_xb       = new TH2F("h2_y_xb"       ,"protons (p/E_{miss} cuts) (Full angular range);y [GeV];x_{B};Counts",200,-0.2,0.35,200,0.4,1.25);

	float l_Q2_min = 1.;
	float l_Q2_max = 4.;

	TH2F * h2_nu_Q2_ec   = new TH2F("h2_nu_Q2_ec"   ,"protons (p/E_{miss} cuts) (EC angular range);#omega [GeV];Q^{2} [GeV^{2}];Counts"
			,  200 ,  0.0 , 2.20 ,  200 , l_Q2_min , l_Q2_max );
	TH2F * h2_nu_Q2_full = new TH2F("h2_nu_Q2_full" ,"protons (p/E_{miss} cuts) (Full angular range);#omega [GeV];Q^{2} [GeV^{2}];Counts"
			,  200 ,  0.0 , 2.20 ,  200 , 1.0 , 4.00 );
	TH1F * h1_xb_ec      = new TH1F("h1_xb_ec"      ,"protons (p/E_{miss} cuts);x_{B};Counts"                       , 200 ,  0.1 , 1.80 );
	TH1F * h1_xb_full    = new TH1F("h1_xb_full"    ,"protons (p/E_{miss} cuts);x_{B};Counts"                       , 200 ,  0.1 , 1.80 );
	h1_xb_ec -> SetLineColor(2);

	// ----------------------------------
	// Checking the QE event selection
	float l_pe_min = 2.9;
	float l_pe_max = 3.6;
	TH2F * h2_pe_ppsm    = new TH2F("h2_pe_ppsm"   ,"(a) smeared protons (acc + p/E^{opt}_{miss} + e kinem cuts);p_{e} [GeV];p^{sm}_{p} [GeV];Counts"
			,40,l_pe_min,l_pe_max, 40, 1.2 , 2.4 );
	TH2F * h2_pe_pn      = new TH2F("h2_pe_pn"     ,"(b) neutrons (acc + p/E^{opt}_{miss} + e kinem cuts);p_{e} [GeV];p_{n} [GeV];Counts"
			,40,l_pe_min,l_pe_max,40,1.2,2.4);
	TH2F * h2_pe_pp      = new TH2F("h2_pe_pp"     ,"(c) protons (acc + p/E_{miss} + e kinem cuts);p_{e} [GeV];p_{p} [GeV];Counts"       
			,40,l_pe_min,l_pe_max,40,1.2,2.4);

	TH1F * h1_phi_epsm   = new TH1F("h1_phi_epsm"  ,"(a) smeared protons (acc + p/E^{opt}_{miss} + e kinem cuts);|#phi^{sm}_{p}-#phi_{e}| [deg];Counts" , 40, 160, 200);
	TH1F * h1_phi_en     = new TH1F("h1_phi_en"    ,"(b) neutrons (acc + p/E^{opt}_{miss} + e kinem cuts);|#phi_{n}-#phi_{e}| [deg];Counts"             , 40, 160, 200);
	TH1F * h1_phi_ep     = new TH1F("h1_phi_ep"    ,"(c) protons (acc + p/E_{miss} + e kinem cuts);|#phi_{p}-#phi_{e}| [deg];Counts"              , 40, 160, 200);

	TH1F * h1_theta_e_ps = new TH1F("h1_theta_e_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);#theta^{measured}_{e}-#theta^{expected}_{e} [deg];Counts",40,-25,30);
	TH1F * h1_theta_e_n  = new TH1F("h1_theta_e_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta^{measured}_{e}-#theta^{expected}_{e} [deg];Counts",40,-25,30);
	TH1F * h1_theta_ps   = new TH1F("h1_theta_ps"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta^{measured}_{N}-#theta^{expected}_{N} [deg];Counts",40,-25,30);
	TH1F * h1_theta_n    = new TH1F("h1_theta_n"   ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta^{measured}_{N}-#theta^{expected}_{N} [deg];Counts",40,-25,30);
	h1_theta_e_ps -> SetLineColor(2);
	h1_theta_ps   -> SetLineColor(2);

	float l_th_pq_min = 0.82;
	float l_th_pq_max = 1.22;
	TH2F * h2_thpq_ppq   = new TH2F("h2_thpq_ppq"  ,"(a) smeared protons (acc + p/E^{opt}_{miss} + e kinem cuts);p_{p}/q;#theta_{pq};Counts"
			,40,l_th_pq_min,l_th_pq_max,40,0,12);
	TH2F * h2_thnq_pnq   = new TH2F("h2_thnq_pnq"  ,"(b) neutrons (acc + p/E^{opt}_{miss} + e kinem cuts);p_{n}/q;#theta_{nq};Counts"
			,40,l_th_pq_min,l_th_pq_max,40,0,12);

	TH1F * h1_p_e_ps  = new TH1F("h1_p_e_ps" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);p_{e} [GeV];Counts"	  , 40 , 2.9 , 3.8 );
	TH1F * h1_p_e_n   = new TH1F("h1_p_e_n"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);p_{e} [GeV];Counts"	  , 40 , 2.9 , 3.8 );
	TH1F * h1_y_e_ps  = new TH1F("h1_y_e_ps" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);y [GeV];Counts"	  , 40 ,-0.1 , 0.3 );
	TH1F * h1_y_e_n   = new TH1F("h1_y_e_n"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);y [GeV];Counts"	  , 40 ,-0.1 , 0.3 );
	TH1F * h1_th_e_ps = new TH1F("h1_th_e_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{e} [deg];Counts" , 40 , 14. , 28. );
	TH1F * h1_th_e_n  = new TH1F("h1_th_e_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{e} [deg];Counts" , 40 , 14. , 28. );
	TH1F * h1_Q2_e_ps = new TH1F("h1_Q2_e_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);Q^{2} [GeV^{2}];Counts"  , 40 ,  0. , 3.5 );
	TH1F * h1_Q2_e_n  = new TH1F("h1_Q2_e_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);Q^{2} [GeV^{2}];Counts"  , 40 ,  0. , 3.5 );
	TH1F * h1_Nu_e_ps = new TH1F("h1_Nu_e_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);#omega [GeV];Counts"	  , 40 , 0.7 , 1.7 );
	TH1F * h1_Nu_e_n  = new TH1F("h1_Nu_e_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#omega [GeV];Counts"	  , 40 , 0.7 , 1.7 );
	TH1F * h1_Xb_e_ps = new TH1F("h1_Xb_e_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);x_{B};Counts"		  , 40 , 0.5 , 1.2 );
	TH1F * h1_Xb_e_n  = new TH1F("h1_Xb_e_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);x_{B};Counts"		  , 40 , 0.5 , 1.2 );
	h1_p_e_ps  -> SetLineColor(2);
	h1_y_e_ps  -> SetLineColor(2);
	h1_th_e_ps -> SetLineColor(2);
	h1_Q2_e_ps -> SetLineColor(2);
	h1_Nu_e_ps -> SetLineColor(2);
	h1_Xb_e_ps -> SetLineColor(2);

	TH1F * h1_p_ps   = new TH1F("h1_p_ps"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);p_{N} [GeV];Counts"              , 40 , 1.2 , 2.4 );
	TH1F * h1_p_n    = new TH1F("h1_p_n"   ,"(acc + p/E^{opt}_{miss} + e kinem cuts);p_{N} [GeV];Counts"              , 40 , 1.2 , 2.4 );
	TH1F * h1_th_ps  = new TH1F("h1_th_ps" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{N} [deg];Counts"         , 40 , 30. , 55. );
	TH1F * h1_th_n   = new TH1F("h1_th_n"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{N} [deg];Counts"         , 40 , 30. , 55. );
	TH1F * h1_pq_ps  = new TH1F("h1_pq_ps" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{Nq} [deg];Counts"        , 40 ,  0. ,  8. );
	TH1F * h1_pq_n   = new TH1F("h1_pq_n"  ,"(acc + p/E^{opt}_{miss} + e kinem cuts);#theta_{Nq} [deg];Counts"        , 40 ,  0. ,  8. );
	TH1F * h1_phi_ps = new TH1F("h1_phi_ps","(acc + p/E^{opt}_{miss} + e kinem cuts);|#phi_{N}-#phi_{q}| [deg];Counts", 40 , 160 , 200 );
	TH1F * h1_phi_n  = new TH1F("h1_phi_n" ,"(acc + p/E^{opt}_{miss} + e kinem cuts);|#phi_{N}-#phi_{q}| [deg];Counts", 40 , 160 , 200 );
	h1_p_ps   -> SetLineColor(2);
	h1_th_ps  -> SetLineColor(2);
	h1_pq_ps  -> SetLineColor(2);
	h1_phi_ps -> SetLineColor(2);

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Testing smeared quantities
	double pm_sm, em_sm, pm, em, ph, ph_sm;
	TFile * output_root = new TFile("./output.root","RECREATE");
	TTree * tree = new TTree("t","false pos neg");
	tree  ->Branch("ph"    , &ph     , "ph/D"   );
	tree  ->Branch("ph_sm" , &ph_sm  , "ph_sm/D");
	tree  ->Branch("pm_sm" , &pm_sm  , "pm_sm/D");
	tree  ->Branch("em_sm" , &em_sm  , "em_sm/D");
	tree  ->Branch("pm"    , &pm     , "pm/D"   );
	tree  ->Branch("em"    , &em     , "em/D"   );

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
		h_mom = P4_Coulomb_Corrected(h_mom,Z,Z+N);	// Applying Coulomb corrections

		q4    = V4_E0-e_mom;
		Pmiss = (h_mom - q4).Rho();			// Missing momentum
		Emiss = fn_Emiss(Pmiss,Nu,Mtar,h_mom.E(),Mp);	// Missing energy

		if     (tab_targ == "4He")	Mres = M3H ;
		else if(tab_targ == "3He")	Mres = M2H ;
		else if(tab_targ == "12C")      Mres = M11B;

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
		if(!p_points_to_EC_fid(h_mom.Vect(),vtx_z[1],10.)) continue;
		h1_ep_EC_x_cut    -> Fill(p_ecx      );
		h1_ep_EC_y_cut    -> Fill(p_ecy      );
		h2_ep_EC_xy_cut   -> Fill(p_ecx,p_ecy);
		// ---------------------------------------
		// phi-theta fiducial cut
		theta = h_mom.Vect().Theta()*180./pi;
		phi   = h_mom.Vect().Phi()  *180./pi;
		if(phi<-30) phi+=360;
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
		// Smearing protons to simulate neutrons
		p_smeared = gRandom -> Gaus(h_mom.Vect().Mag(),f_delta_p(h_mom.Vect().Mag()));
		if(p_smeared>2.34) continue;
		pp_min_psm = h_mom.Vect().Mag() - p_smeared;
		h1_pp_mi_psm -> Fill(pp_min_psm);

		// Copying direction of p momentum vector
		u1 = h_mom.Vect().Unit();
		p_sm_mom.SetXYZM(p_smeared*u1.X(),p_smeared*u1.Y(),p_smeared*u1.Z(),Mp);

		Pmiss_sm = (p_sm_mom - q4).Rho();			// Missing momentum
		Emiss_sm = fn_Emiss(Pmiss_sm,Nu,Mtar,p_sm_mom.E(),Mp);  // Missing energy	

		h2_p_Pm_Em   -> Fill(Pmiss,Emiss);
		h2_psm_Pm_Em -> Fill(Pmiss_sm,Emiss_sm);
		h2_p_pm_em   -> Fill(Pmiss   ,Emiss   );
		h2_psm_pm_em -> Fill(Pmiss_sm,Emiss_sm);
		// ------------------------------------------------------------------------------
		// Quasi Elastic cuts	
		y_sc = y_scale_func(q4,Mtar,Mres,Mp);
		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;

		h2_p_Nu_y    -> Fill(Nu,y_sc    );
		h2_p_Q2_thpq -> Fill(Q2,theta_pq);	

		h2_but_pm_em_0 -> Fill(Emiss_sm-Emiss,Pmiss_sm-Pmiss);

		if((Pmiss<Pmiss_cut)&&(Emiss<0.15)&&(Emiss>Emiss_cut)){
                        h2_p_Nu_y_QEcut2    -> Fill(Nu,y_sc    );
                        h2_p_Q2_thpq_QEcut2 -> Fill(Q2,theta_pq);
		}

		if((Pmiss<Pmiss_cut)&&(Emiss<Emiss_cut)){
			h2_p_Nu_y_QEcut    -> Fill(Nu,y_sc    );
			h2_p_Q2_thpq_QEcut -> Fill(Q2,theta_pq);

			if(		(y_sc     > y_scl_min)&&(y_sc < y_scl_max)&&
					(Nu       > Nu_min   )&&(Nu   < Nu_max   )&&
					//(Q2       > Q2_min   )&&(Q2   < Q2_max   )&&
					(theta_pq < Th_pq_max)
			  ){
				h1_y_cut2     -> Fill( y_sc     );
				h1_Nu_cut2    -> Fill( Nu       );
				h1_th_pq_cut2 -> Fill( theta_pq );
				h1_Q2_cut2    -> Fill( Q2       );

				h2_psm_Pm_Em_cut2 -> Fill(Pmiss_sm,Emiss_sm);
			}
		}
		// ----	
		if(		(y_sc     > y_scl_min)&&(y_sc < y_scl_max)&&
				(Nu       > Nu_min   )&&(Nu   < Nu_max   )&&
				//(Q2       > Q2_min   )&&(Q2   < Q2_max   )&&
				(theta_pq < Th_pq_max)
		  ){
			h2_p_Pm_Em_cut1   -> Fill(Pmiss   ,Emiss   );
			h2_psm_Pm_Em_cut1 -> Fill(Pmiss_sm,Emiss_sm);

			h2_p_pm_em_1   -> Fill(Pmiss   ,Emiss   );
			h2_psm_pm_em_1 -> Fill(Pmiss_sm,Emiss_sm);

			h2_pm_pmsm   -> Fill(Pmiss             ,Pmiss_sm );
			h2_em_emsm   -> Fill(Emiss             ,Emiss_sm );
			h2_ph_phsm   -> Fill(h_mom.Vect().Mag(),p_smeared);
			h2_ph_pmsm   -> Fill(h_mom.Vect().Mag(),Pmiss_sm );

			h2_pm_pmsm_pm_diff -> Fill(Pmiss_sm-Pmiss              ,Pmiss             );
			h2_em_emsm_em_diff -> Fill(Emiss_sm-Emiss              ,Emiss             );
			h2_ph_phsm_ph_diff -> Fill(p_smeared-h_mom.Vect().Mag(),h_mom.Vect().Mag());
			h2_ph_pmsm_ph_diff -> Fill(Pmiss_sm-h_mom.Vect().Mag() ,h_mom.Vect().Mag());

			h2_but_pm_em_1 -> Fill(Emiss_sm-Emiss              ,Pmiss_sm-Pmiss);
			h2_but_pm_ph   -> Fill(p_smeared-h_mom.Vect().Mag(),Pmiss_sm-Pmiss);
			h2_but_em_ph   -> Fill(p_smeared-h_mom.Vect().Mag(),Emiss_sm-Emiss);

			// Filling output tree
			ph    = h_mom.Vect().Mag();
			ph_sm = p_sm_mom.Vect().Mag();
			pm_sm = Pmiss_sm;
			em_sm = Emiss_sm;
			pm    = Pmiss   ;
			em    = Emiss   ;
			tree -> Fill()  ;

		}

		// False positive and False negative rates
		if(             (y_sc     > y_scl_min)&&(y_sc < y_scl_max)&&
				(Nu       > Nu_min   )&&(Nu   < Nu_max   )&&
				//(Q2       > Q2_min   )&&(Q2   < Q2_max   )&&
				(theta_pq < Th_pq_max)
		  ){
			for(int i = 0 ; i < no_pmiss_cuts ; i++){
				for(int j = 0 ; j < no_emiss_cuts ; j++){
					// False positive
					if((Pmiss_sm<Pmiss_cut_lim[i])&&(Emiss_sm<Emiss_cut_lim[j])){
						fal_denom[j][i]++;
						if(!((Pmiss<Pmiss_cut)&&(Emiss<Emiss_cut))){
							fal_pos_num[j][i]++;
						}
					}
					// -------
					// False negative
					if((Pmiss<Pmiss_cut)&&(Emiss<Emiss_cut)){
						if(!((Pmiss_sm<Pmiss_cut_lim[i])&&(Emiss_sm<Emiss_cut_lim[j]))){ 
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

	output_root->cd   ();
	tree       ->Write();
	output_root->Close();

	double A, B, eA, eB;

	TGraphErrors ** f_fal_pos = new TGraphErrors*[no_emiss_cuts];
	TLegend * leg_false_pos = new TLegend(0.20,0.50,0.40,0.85);
	leg_false_pos -> SetLineColor(0);
	leg_false_pos -> SetHeader("E_{miss} cut [GeV]");
	for(int j = 0 ; j < no_emiss_cuts ; j++){
		for(int i = 0 ; i < no_pmiss_cuts ; i++){
			A  = fal_pos_num[j][i];
			B  = fal_denom  [j][i];
			eA = sqrt(A); 
			eB = sqrt(B);
			fal_pos  [j][i] = 100.*A/B;
			fal_pos_e[j][i] = 100.*sqrt(eA*eA - A*A*eB*eB/(B*B))/B;
		}
		f_fal_pos[j] = new TGraphErrors(no_pmiss_cuts,Pmiss_cut_lim,fal_pos[j],0,fal_pos_e[j]);
		f_fal_pos[j] -> SetMarkerColor(j+1);
		f_fal_pos[j] -> SetMarkerSize (1  );
		f_fal_pos[j] -> SetMarkerStyle(20 );
		f_fal_pos[j] -> SetFillColor  (0  );
		f_fal_pos[j] -> SetMinimum(  0);
		f_fal_pos[j] -> SetMaximum(70);
		f_fal_pos[j] -> SetTitle("");
		f_fal_pos[j] -> GetXaxis() -> SetTitle("p_{miss} cut [GeV]");
		f_fal_pos[j] -> GetYaxis() -> SetTitle("False Positive [%]");
		leg_false_pos -> AddEntry(f_fal_pos[j],Form("%.2f",Emiss_cut_lim[j]));
	}
	// False negative rate
	TGraphErrors ** f_fal_neg = new TGraphErrors*[no_emiss_cuts];
	TLegend * leg_false_neg = new TLegend(0.20,0.20,0.40,0.55);
	leg_false_neg -> SetLineColor(0);
	leg_false_neg -> SetHeader("E_{miss} cut [GeV]");
	for(int j = 0 ; j < no_emiss_cuts ; j++){
		for(int i = 0 ; i < no_pmiss_cuts ; i++){
			A  = fal_neg_num[j][i];	
			B  = fal_denom  [j][i];
			eA = sqrt(A);
			eB = sqrt(B);
			fal_neg  [j][i] = 100.*A/B;
			fal_neg_e[j][i] = 100.*sqrt(eA*eA - A*A*eB*eB/(B*B))/B;
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
		leg_false_neg -> AddEntry(f_fal_neg[j],Form("%.2f",Emiss_cut_lim[j]));
	}
	// False negative - false positive
	TGraphErrors ** f_pos_neg = new TGraphErrors*[no_emiss_cuts];
	TLegend * leg_pos_neg = new TLegend(0.65,0.50,0.85,0.85);
	leg_pos_neg -> SetLineColor(0);
	leg_pos_neg -> SetHeader("E_{miss} cut [GeV]");
	for(int j = 0 ; j < no_emiss_cuts ; j++){ 
		f_pos_neg[j] = new TGraphErrors(no_pmiss_cuts,fal_pos[j],fal_neg[j],fal_pos_e[j],fal_neg_e[j]);
		f_pos_neg[j] -> SetMarkerColor(j+1);
		f_pos_neg[j] -> SetMarkerSize (1  );
		f_pos_neg[j] -> SetMarkerStyle(20 );
		f_pos_neg[j] -> SetFillColor  (0  );
		f_pos_neg[j] -> SetMinimum(  0);
		f_pos_neg[j] -> SetMaximum(100);
		f_pos_neg[j] -> GetXaxis() -> SetLimits(0,100);
		f_pos_neg[j] -> SetTitle("");
		f_pos_neg[j] -> GetXaxis() -> SetTitle("False Positive [%]");
		f_pos_neg[j] -> GetYaxis() -> SetTitle("False Negative [%]");
		leg_pos_neg -> AddEntry(f_pos_neg[j],Form("%.2f",Emiss_cut_lim[j]));
	}
	// Optimizing cut
	int min_i     = 99999;
	int min_j     = 99999;
	float min_val = 99999;
	float tem_val = 99999;
	for(int i = 0 ; i < no_pmiss_cuts ; i++){
		for(int j = 0 ; j < no_emiss_cuts ; j++){
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
	min_i = 10;
	min_j =  3;

	cout << "**************************\nThe simultaneous minimization of false positive and negative rate happens for:\n" << 
		"Pmiss_cut = " << Pmiss_cut_lim[min_i] << " GeV, and\nEmiss_cut = " << Emiss_cut_lim[min_j] << " GeV\n" <<
		"and corresponds to a false positive of: " << fal_pos[min_j][min_i] << "%\nand a false negative of: "<<
		fal_neg[min_j][min_i] << "%\n";

	// Creating circles to point to where the optimized cuts are
	TEllipse *el1 = new TEllipse(Pmiss_cut_lim[min_i] ,fal_pos[min_j][min_i],.008,1.0);
	el1->SetFillColor(0);
	el1->SetFillStyle(0);

	TEllipse *el2 = new TEllipse(Pmiss_cut_lim[min_i] ,fal_neg[min_j][min_i],.008,3.0);
	el2->SetFillColor(0);
	el2->SetFillStyle(0);

	TEllipse *el3 = new TEllipse(fal_pos[min_j][min_i],fal_neg[min_j][min_i], 0.8,3.0);
	el3->SetFillColor(0);
	el3->SetFillStyle(0);

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

		if     (tab_targ == "4He")      Mres = M3H ;
		else if(tab_targ == "3He")      Mres = M2H ;
		else if(tab_targ == "12C")      Mres = M11B;

		y_sc = y_scale_func(q4,Mtar,Mres,Mp);  
		Xb   = Q2/(2.*Nu*Mp);

		theta = h_mom.Vect().Theta()*180./pi;
		phi   = h_mom.Vect().Phi()  *180./pi;
		if(phi<-30) phi+=360;

		theta_e = e_mom.Vect().Theta()*180./pi;

		phi_e = e_mom.Vect().Phi()  *180./pi;
		if(phi_e<-30) phi_e+=360;

		Pmiss_par = (h_mom - q4).Vect()*q4.Vect().Unit();

		theta_pq = acos(h_mom.Vect().Unit()*q4.Vect().Unit())*180./pi;
		// ---------------------------------------
		// Momentum cut
		if(h_mom.Vect().Mag()>2.34) continue;

		// ------------------------------------------------------------------------------
		// Quasi Elastic cuts
		if((Pmiss<Pmiss_cut)&&(Emiss<Emiss_cut)){  
			h2_y_xb       -> Fill(y_sc ,Xb       );
			h2_thp_y      -> Fill(theta,y_sc     );
			h1_y_full     -> Fill(y_sc           );
			h2_thp_pmpar  -> Fill(theta,Pmiss_par);
			h1_pmpar_full -> Fill(Pmiss_par      );
			h1_xb_full    -> Fill(Xb             );
			h2_nu_Q2_full -> Fill(Nu   ,Q2       );	
			if(theta<44.){	// Cutting top of EC
				h1_y_ec     -> Fill(y_sc     );
				h1_pmpar_ec -> Fill(Pmiss_par);
				h1_xb_ec    -> Fill(Xb       );
				h2_nu_Q2_ec -> Fill(Nu,Q2    );
			}
		}

		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		if(!p_points_to_EC_fid(h_mom.Vect(),vtx_z[1],10.)) continue;	// EC 10 cm cut 
		if(theta>44.) continue;						// Cutting top of EC
		// Momentum cut - Already applied this cut above

		// ------------------------------------------------------------------------------
		// QE cuts
		if((Pmiss<Pmiss_cut)&&(Emiss<Emiss_cut)){
			if(             (y_sc     > y_scl_min)&&(y_sc     < y_scl_max)&&
					(Nu       > Nu_min   )&&(Nu       < Nu_max   )&&
					(Q2       > Q2_min   )&&(Q2       < Q2_max   )&&
					(theta_pq < Th_pq_max)  
			  ){
				h2_pe_pp   -> Fill(e_mom.Rho(),h_mom   .Rho());	
				h1_phi_ep  -> Fill(phi-phi_e);

				if(p_points_to_EC_fid((q4+V4_Targ).Vect(),vtx_z[0],10.)){
					nu_eep      += 1.;
					raw_nu_p    += 1.;
					nu_eep_een_GK05 += (sigma_p_GK05->Eval(Q2,theta_e))/(sigma_n_GK05->Eval(Q2,theta_e));
					nu_eep_een_AMT  += (sigma_p_AMT ->Eval(Q2,theta_e))/(sigma_n_AMT ->Eval(Q2,theta_e));
					nu_eep_een_MMD  += (sigma_p_MMD ->Eval(Q2,theta_e))/(sigma_n_MMD ->Eval(Q2,theta_e));

					for(int i1 = 0 ; i1 < bins_eepn ; i1++){
						if(h_mom.Rho()>p_bins_eepn[i1]&&h_mom.Rho()<p_bins_eepn[i1+1])
							nu_eep_p[i1] += 1.;	
					}
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

		if     (tab_targ == "4He")      Mres = M3H ;
		else if(tab_targ == "3He")      Mres = M2H ;
		else if(tab_targ == "12C")      Mres = M11B;

		// ------------------------------------------------------------------------------ 
		// Smearing protons to simulate neutrons
		p_smeared = gRandom -> Gaus(h_mom.Vect().Mag(),f_delta_p(h_mom.Vect().Mag()));
		if(p_smeared>2.34) continue;
		pp_min_psm = h_mom.Vect().Mag() - p_smeared;

		// Copying direction of p momentum vector
		u1 = h_mom.Vect().Unit();
		p_sm_mom.SetXYZM(p_smeared*u1.X(),p_smeared*u1.Y(),p_smeared*u1.Z(),Mp);

		Pmiss_sm = (p_sm_mom - q4).Rho();                       // Missing momentum
		Emiss_sm = fn_Emiss(Pmiss_sm,Nu,Mtar,p_sm_mom.E(),Mp);  // Missing energy
		// -----------------------
		theta = p_sm_mom.Vect().Theta()*180./pi;
		phi   = p_sm_mom.Vect().Phi()  *180./pi;
		if(phi<-30) phi+=360;

		theta_e = e_mom.Vect().Theta()*180./pi;
		phi_e   = e_mom.Vect().Phi()  *180./pi;
		if(phi_e<-30) phi_e+=360;

		theta_e_calc = acos((Ebeam - p_sm_mom.Rho()*cos(theta  *pi/180.))/e_mom.   Rho())*180./pi;
		theta_h_calc = acos((Ebeam - e_mom.   Rho()*cos(theta_e*pi/180.))/p_sm_mom.Rho())*180./pi;
		// ---------------------------------------
		// Momentum cut
		if(p_sm_mom.Vect().Mag()>2.34) continue;
		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		if(!p_points_to_EC_fid(p_sm_mom.Vect(),vtx_z[1],10.)) continue;	// EC 10 cm cut
		if(theta>44.) continue;						// Cutting top of EC
		// Momentum cut - Already applied this cut above

		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;
		p_over_q = p_sm_mom.Rho()/q4.Rho();
		y_sc = y_scale_func(q4,Mtar,Mres,Mp);
		// ------------------------------------------------------------------------------
		// QE cuts 
		if((Pmiss_sm<Pmiss_cut_lim[min_i])&&(Emiss_sm<Emiss_cut_lim[min_j])){
			if(             (y_sc     > y_scl_min)&&(y_sc     < y_scl_max)&&
					(Nu       > Nu_min   )&&(Nu       < Nu_max   )&&
					(Q2       > Q2_min   )&&(Q2       < Q2_max   )
			  ){
				h2_thpq_ppq   -> Fill(p_over_q,theta_pq);

				if (theta_pq < Th_pq_max) {
					h2_pe_ppsm    -> Fill(e_mom.Rho(),p_sm_mom.Rho());
					h1_phi_epsm   -> Fill(phi-phi_e);
					h1_theta_e_ps -> Fill(theta_e-theta_e_calc);
					h1_theta_ps   -> Fill(theta  -theta_h_calc);

					h1_p_e_ps  -> Fill(e_mom.Rho()  );
					h1_y_e_ps  -> Fill(y_sc         );
					h1_th_e_ps -> Fill(theta_e      );
					h1_Q2_e_ps -> Fill(Q2           );
					h1_Nu_e_ps -> Fill(Nu           );
					h1_Xb_e_ps -> Fill(Q2/(2.*Nu*Mp));

					h1_p_ps   -> Fill(p_sm_mom.Rho());
					h1_th_ps  -> Fill(theta         );
					h1_pq_ps  -> Fill(theta_pq      );
					h1_phi_ps -> Fill(phi-phi_e     );
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

		if     (tab_targ == "4He")      Mres = M3H ;
		else if(tab_targ == "3He")      Mres = 2*Mp;
		else if(tab_targ == "12C")      Mres = M11C;

		// ---------------------------------------
		// Weight from neutron detection eff.
		if(h_mom.Vect().Mag()<0.5) continue;
		n_weight = 1./n_efficiency(h_mom.Rho());

		// ------------------------------------------------------------------------------
		// Matching proton and neutron acceptances
		// ---------------------------------------
		// phi-theta fiducial cut
		theta = h_mom.Vect().Theta()*180./pi;
		phi   = h_mom.Vect().Phi()  *180./pi;
		if(phi<-30) phi+=360;

		theta_e = e_mom.Vect().Theta()*180./pi;
		phi_e   = e_mom.Vect().Phi()  *180./pi;
		if(phi_e<-30) phi_e+=360;

		theta_e_calc = acos((Ebeam - h_mom.Rho()*cos(theta  *pi/180.))/e_mom.Rho())*180./pi;
		theta_h_calc = acos((Ebeam - e_mom.Rho()*cos(theta_e*pi/180.))/h_mom.Rho())*180./pi;

		h1_en_th     -> Fill(    theta);
		h1_en_phi    -> Fill(      phi);
		h2_en_th_phi -> Fill(phi,theta);

		if(theta>44.) continue;	// Cutting top of EC
		if(!pFiducialCut(h_mom.Vect())) continue;
		h1_en_th_cut     -> Fill(    theta );
		h1_en_phi_cut    -> Fill(      phi );
		h2_en_th_phi_cut -> Fill(phi,theta );

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
		h2_n_Pm_Em -> Fill(Pmiss  , Emiss , n_weight);
		h2_n_pm_em -> Fill(Pmiss  , Emiss , n_weight);

		// ----
		// Copying direction of p momentum vector
		u1 = h_mom.Vect().Unit();

		y_sc = y_scale_func(q4,Mtar,Mres,Mn);
		theta_pq = acos(u1*q4.Vect().Unit())*180./pi;
		p_over_q = h_mom.Rho()/q4.Rho();

		if(		(y_sc     > y_scl_min)&&(y_sc     < y_scl_max)&&
				(Nu       > Nu_min   )&&(Nu       < Nu_max   )&&
				(Q2       > Q2_min   )&&(Q2       < Q2_max   )&&
				(theta_pq < Th_pq_max)
		  ){
			h2_n_Pm_Em_cut1   -> Fill(Pmiss   ,Emiss   , n_weight);
			h2_n_pm_em_1      -> Fill(Pmiss   ,Emiss   , n_weight);
		}

		// ------------------------------------------------------------------------------
		// QE cuts
		if((Pmiss<Pmiss_cut_lim[min_i])&&(Emiss<Emiss_cut_lim[min_j])){
			if(             (y_sc     > y_scl_min)&&(y_sc     < y_scl_max)&&
					(Nu       > Nu_min   )&&(Nu       < Nu_max   )&&
					(Q2       > Q2_min   )&&(Q2       < Q2_max   ) 
			  ){
				h2_thnq_pnq  -> Fill(p_over_q,theta_pq, n_weight);        

				if (theta_pq < Th_pq_max) {
					h2_pe_pn     -> Fill(e_mom.Rho(),h_mom.Rho(), n_weight);
					h1_phi_en    -> Fill(phi-phi_e              , n_weight);
					h1_theta_e_n -> Fill(theta_e-theta_e_calc   , n_weight);
					h1_theta_n   -> Fill(theta  -theta_h_calc   , n_weight);

					h1_p_e_n  -> Fill(e_mom.Rho()  , n_weight);
					h1_y_e_n  -> Fill(y_sc         , n_weight);
					h1_th_e_n -> Fill(theta_e      , n_weight);
					h1_Q2_e_n -> Fill(Q2           , n_weight);
					h1_Nu_e_n -> Fill(Nu           , n_weight);
					h1_Xb_e_n -> Fill(Q2/(2.*Nu*Mn), n_weight);

					h1_p_n    -> Fill(h_mom.Rho()  , n_weight);
					h1_th_n   -> Fill(theta        , n_weight);
					h1_pq_n   -> Fill(theta_pq     , n_weight);
					h1_phi_n  -> Fill(phi-phi_e    , n_weight);

					if(p_points_to_EC_fid((q4+V4_Targ).Vect(),vtx_z[0],10.)){

						nu_een      += n_weight;
						raw_nu_n    += 1.;

						for(int i1 = 0 ; i1 < bins_eepn ; i1++){
							if(h_mom.Rho()>p_bins_eepn[i1]&&h_mom.Rho()<p_bins_eepn[i1+1])
								nu_een_p[i1] += n_weight;
						}
					}
				}
			}
		}

	}
	cout << "**************************\n\n";

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Preparing plots
	TLegend ** leg_ec_cut = new TLegend*[4];
	for(int i = 0 ; i < 4 ; i++){
		leg_ec_cut[i] = new TLegend(0.2,0.7,0.4,0.85);
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

	h2_ep_EC_xy -> SetTitle("protons before and after 10 cm EC cut"  );
	h2_en_EC_xy -> SetTitle("protons and neutrons after 10 cm EC cut");

	TLegend ** leg_th_phi_cut = new TLegend*[4];
	for(int i = 0 ; i < 4 ; i++){
		leg_th_phi_cut[i] = new TLegend(0.2,0.7,0.4,0.85);
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
	TLegend * leg_p_cut = new TLegend(0.2,0.7,0.4,0.85);
	leg_p_cut -> SetLineColor(0);
	leg_p_cut -> AddEntry(h1_p_p    ,"p evts before cut");
	leg_p_cut -> AddEntry(h1_p_p_cut,"p evts after  cut");
	leg_p_cut -> AddEntry(h1_n_p    ,"n evts before cut");
	leg_p_cut -> AddEntry(h1_n_p_cut,"n evts after  cut");
	// ----
	//TBox *b_QE_cuts = new TBox(0.0,0.0,Pmiss_cut,Emiss_cut);
	//b_QE_cuts -> SetFillStyle(0);
	//b_QE_cuts -> SetFillColor(0);
	//b_QE_cuts -> SetLineColor(2);
	//b_QE_cuts -> SetLineWidth(3);
	float low_em_bound = h2_p_Pm_Em -> GetYaxis() -> GetBinLowEdge(1);
	TLine * l1_QE_cuts = new TLine(0.0      ,Emiss_cut   ,Pmiss_cut,Emiss_cut);
	TLine * l2_QE_cuts = new TLine(Pmiss_cut,low_em_bound,Pmiss_cut,Emiss_cut);
	l1_QE_cuts -> SetLineColor(2);
	l1_QE_cuts -> SetLineWidth(3);
	l2_QE_cuts -> SetLineColor(2);
	l2_QE_cuts -> SetLineWidth(3);
	// ----
	TBox *b_Nu_y = new TBox(Nu_min,y_scl_min,Nu_max,y_scl_max);
	b_Nu_y -> SetFillStyle(0);
	b_Nu_y -> SetFillColor(0);
	b_Nu_y -> SetLineColor(2);
	b_Nu_y -> SetLineWidth(3);
	// ----
	TBox *b_Q2_th_pq = new TBox(Q2_min,0.0,Q2_max,Th_pq_max);
	b_Q2_th_pq -> SetFillStyle(0);
	b_Q2_th_pq -> SetFillColor(0);
	b_Q2_th_pq -> SetLineColor(2);
	b_Q2_th_pq -> SetLineWidth(3);
	// ----
	TLegend * leg_ydist_eccut = new TLegend(0.25,0.7,0.40,0.85);
	leg_ydist_eccut -> SetLineColor(0);
	leg_ydist_eccut -> AddEntry(h1_y_full, "Without EC ang region cut");
	leg_ydist_eccut -> AddEntry(h1_y_ec  , "With EC ang region cut"   );
	// ----
	TLine * l_ec_ang_reg1 = new TLine(45,-0.5,45,0.6);
	l_ec_ang_reg1 -> SetLineColor(2);
	l_ec_ang_reg1 -> SetLineWidth(2);
	// ----
	TLine * l_ec_ang_reg2 = new TLine(45,-0.3,45,0.3);
	l_ec_ang_reg2 -> SetLineColor(2);
	l_ec_ang_reg2 -> SetLineWidth(2);
	// ----
	float slp = 2*Mp;
	float l_Nu_min = l_Q2_min/slp;
	float l_Nu_max = l_Q2_max/slp;
	TLine * l_Q2_Nu = new TLine(l_Nu_min,l_Q2_min,l_Nu_max,l_Q2_max);
	l_Q2_Nu -> SetLineColor(2);
	l_Q2_Nu -> SetLineWidth(2); 
	// ----
	float l_p_p_min = sqrt( pow( Ebeam+Mp - l_pe_min ,2)-Mp*Mp);
	float l_p_p_max = sqrt( pow( Ebeam+Mp - l_pe_max ,2)-Mp*Mp);
	float l_p_n_min = sqrt( pow( Ebeam+Mn - l_pe_min ,2)-Mn*Mn);
	float l_p_n_max = sqrt( pow( Ebeam+Mn - l_pe_max ,2)-Mn*Mn);

	TLine * l_pe_p = new TLine(l_pe_min,l_p_p_min,l_pe_max,l_p_p_max);
	l_pe_p -> SetLineColor(2);
	l_pe_p -> SetLineWidth(2);
	TLine * l_pe_n = new TLine(l_pe_min,l_p_n_min,l_pe_max,l_p_n_max);
	l_pe_n -> SetLineColor(2);
	l_pe_n -> SetLineWidth(2);

	double temp_integral;
	temp_integral = (h1_theta_e_ps->Integral())/(h1_theta_e_n->Integral());	h1_theta_e_n -> Scale(temp_integral);
	temp_integral = (h1_theta_ps->Integral  ())/(h1_theta_n->Integral  ());	h1_theta_n   -> Scale(temp_integral);

	TLegend * leg_theta_meas_calc = new TLegend(0.7,0.7,0.85,0.85);
	leg_theta_meas_calc -> SetLineColor(0);
	leg_theta_meas_calc -> AddEntry(h1_theta_ps,"Smeared Protons");
	leg_theta_meas_calc -> AddEntry(h1_theta_n ,"Neutrons"       );

	TLine * l_th_pq = new TLine(l_th_pq_min, Th_pq_max, l_th_pq_max, Th_pq_max);
	l_th_pq -> SetLineColor(2);
	l_th_pq -> SetLineWidth(2);

	temp_integral = (h1_p_e_ps->Integral ())/(h1_p_e_n->Integral ());	h1_p_e_n  -> Scale(temp_integral);
	temp_integral = (h1_th_e_ps->Integral())/(h1_th_e_n->Integral());	h1_th_e_n -> Scale(temp_integral);
	temp_integral = (h1_Q2_e_ps->Integral())/(h1_Q2_e_n->Integral());	h1_Q2_e_n -> Scale(temp_integral);
	temp_integral = (h1_Nu_e_ps->Integral())/(h1_Nu_e_n->Integral());	h1_Nu_e_n -> Scale(temp_integral);
	temp_integral = (h1_y_e_ps->Integral ())/(h1_y_e_n->Integral ());	h1_y_e_n  -> Scale(temp_integral);
	temp_integral = (h1_Xb_e_ps->Integral())/(h1_Xb_e_n->Integral());	h1_Xb_e_n -> Scale(temp_integral);

	temp_integral = (h1_p_ps->Integral  ())/(h1_p_n->Integral  ());		h1_p_n   -> Scale(temp_integral);
	temp_integral = (h1_th_ps->Integral ())/(h1_th_n->Integral ());		h1_th_n  -> Scale(temp_integral);
	temp_integral = (h1_pq_ps->Integral ())/(h1_pq_n->Integral ());		h1_pq_n  -> Scale(temp_integral);
	temp_integral = (h1_phi_ps->Integral())/(h1_phi_n->Integral());		h1_phi_n -> Scale(temp_integral);

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Editting plots for presentation
	pretty_TH2F(h2_p_Pm_Em         );
	pretty_TH2F(h2_n_Pm_Em         );
	pretty_TH2F(h2_psm_Pm_Em       );
	pretty_TH2F(h2_p_Nu_y          );
	pretty_TH2F(h2_p_Nu_y_QEcut    );
	pretty_TH2F(h2_p_Nu_y_QEcut2   );
	pretty_TH2F(h2_p_Q2_thpq       );
	pretty_TH2F(h2_p_Q2_thpq_QEcut );
	pretty_TH2F(h2_p_Q2_thpq_QEcut2);
	pretty_TH2F(h2_psm_Pm_Em_cut1  );
	pretty_TH2F(h2_n_Pm_Em_cut1    );
	pretty_TH2F(h2_p_Pm_Em_cut1    );
	pretty_TH2F(h2_psm_Pm_Em_cut2  );
	pretty_TH2F(h2_thp_y           );
	pretty_TH2F(h2_thp_pmpar       );
	pretty_TH2F(h2_y_xb            );
	pretty_TH2F(h2_nu_Q2_full      );
	pretty_TH2F(h2_nu_Q2_ec        );
	pretty_TH2F(h2_pe_ppsm         );
	pretty_TH2F(h2_pe_pn           );
	pretty_TH2F(h2_pe_pp           );
	pretty_TH2F(h2_thpq_ppq        );
	pretty_TH2F(h2_thnq_pnq        );
	pretty_TH2F(h2_but_pm_em_0     );
	pretty_TH2F(h2_but_pm_em_1     );
	pretty_TH2F(h2_but_pm_ph       );
	pretty_TH2F(h2_but_em_ph       );

	pretty_TH1F(h1_pp_mi_psm      );
	pretty_TH1F(h1_y_full         );
	pretty_TH1F(h1_pmpar_full     );
	pretty_TH1F(h1_xb_full        );
	pretty_TH1F(h1_phi_epsm       );
	pretty_TH1F(h1_phi_en         );
	pretty_TH1F(h1_phi_ep         );
	pretty_TH1F(h1_theta_e_ps     );
	pretty_TH1F(h1_theta_ps       );
	pretty_TH1F(h1_p_e_ps         );
	pretty_TH1F(h1_th_e_ps        );
	pretty_TH1F(h1_Q2_e_ps        );
	pretty_TH1F(h1_Nu_e_ps        );
	pretty_TH1F(h1_y_e_ps         );
	pretty_TH1F(h1_Xb_e_ps        );
	pretty_TH1F(h1_p_ps           );
	pretty_TH1F(h1_th_ps          );
	pretty_TH1F(h1_pq_ps          );
	pretty_TH1F(h1_phi_ps         );

	for(int j = 0 ; j < no_emiss_cuts ; j++){
		pretty_TGraphErrors(f_fal_pos[j]);
		pretty_TGraphErrors(f_fal_neg[j]);
		pretty_TGraphErrors(f_pos_neg[j]);
	}

	// -------------------------------------------------------------------------------------------------------------------------------------
	// Plotting histograms
	gStyle->SetOptStat(0);

	TCanvas * c1 = new TCanvas("c1","",1300,800);
	c1 -> Divide(2,2);
	c1 -> cd(1);    h1_en_EC_x  -> Draw();  h1_ep_EC_x_cut  -> Draw("same");        h1_ep_EC_x -> Draw("same");     leg_ec_cut[0]-> Draw("same");
	c1 -> cd(2);    h1_en_EC_y  -> Draw();  h1_ep_EC_y_cut  -> Draw("same");        h1_ep_EC_y -> Draw("same");     leg_ec_cut[1]-> Draw("same");
	c1 -> cd(3);    h2_ep_EC_xy -> Draw();  h2_ep_EC_xy_cut -> Draw("same");                                        leg_ec_cut[2]-> Draw("same");
	c1 -> cd(4);    h2_en_EC_xy -> Draw();  h2_ep_EC_xy_cut -> Draw("same");                                        leg_ec_cut[3]-> Draw("same");
	c1->Modified(); c1->Update();

	TCanvas * c2 = new TCanvas("c2","",1300,800);
	c2 -> Divide(2,2);
	c2 -> cd(1);    h1_en_th     -> Draw(); h1_ep_th         -> Draw("same");       h1_en_th_cut  -> Draw("same");  leg_th_phi_cut[0]-> Draw("same");
	c2 -> cd(2);    h1_en_phi    -> Draw(); h1_ep_phi        -> Draw("same");       h1_en_phi_cut -> Draw("same");  leg_th_phi_cut[1]-> Draw("same");
	c2 -> cd(3);    h2_en_th_phi -> Draw(); h2_en_th_phi_cut -> Draw("same");                                       leg_th_phi_cut[2]-> Draw("same");
	c2 -> cd(4);    h2_ep_th_phi -> Draw(); h2_en_th_phi_cut -> Draw("same");                                       leg_th_phi_cut[3]-> Draw("same");
	c2->Modified(); c2->Update();

	TCanvas * c4 = new TCanvas("c4","",1300,800);
	c4 -> Divide(2,2);
	c4 -> cd(1);    h2_p_Pm_Em  -> Draw("COLZ" );	l1_QE_cuts -> Draw("same");	l2_QE_cuts -> Draw("same");	//b_QE_cuts -> Draw("same");
	c4 -> cd(2);    h2_n_Pm_Em  -> Draw("COLZ" );	l1_QE_cuts -> Draw("same");	l2_QE_cuts -> Draw("same");
	c4 -> cd(3);    h1_pp_mi_psm-> Draw("same" );
	c4 -> cd(4);    h2_psm_Pm_Em-> Draw("COLZ" );
	c4->Modified(); c4->Update();

	TCanvas * c5 = new TCanvas("c5","",1300,800);
	c5 -> Divide(2,2);
	c5 -> cd(1);    h2_p_Nu_y          -> Draw("COL");
	c5 -> cd(2);    h2_p_Nu_y_QEcut    -> Draw("COL");      b_Nu_y -> Draw("same");
	c5 -> cd(3);    h2_p_Q2_thpq       -> Draw("COL");
	c5 -> cd(4);    h2_p_Q2_thpq_QEcut -> Draw("COL");      b_Q2_th_pq -> Draw("same");
	c5->Modified(); c5->Update();

	TCanvas * c6 = new TCanvas("c6","",1300,800);
	c6 -> Divide(2,2);
	c6 -> cd(1);    h1_y_cut2     -> Draw();
	c6 -> cd(2);    h1_Nu_cut2    -> Draw();
	c6 -> cd(3);    h1_th_pq_cut2 -> Draw();
	c6 -> cd(4);    h1_Q2_cut2    -> Draw();
	c6->Modified(); c6->Update();

	TCanvas * c7 = new TCanvas("c7","",1300,800);
	c7 -> Divide(2,2);
	c7 -> cd(1);    h2_psm_Pm_Em_cut1 -> Draw("COLZ");
	c7 -> cd(2);    h2_n_Pm_Em_cut1   -> Draw("COLZ");
	c7 -> cd(3);    h2_p_Pm_Em_cut1   -> Draw("COLZ");
	c7 -> cd(4);    h2_psm_Pm_Em_cut2 -> Draw("COLZ");
	c7->Modified(); c7->Update();

	TCanvas * c8 = new TCanvas("c8","",1300,800);
	c8 -> Divide(2,2);
	c8 -> cd(1);
	f_fal_pos[no_emiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_emiss_cuts-1 ; j++) f_fal_pos[j] -> Draw("samePE");
	leg_false_pos -> Draw("same");
	el1 -> Draw("same");
	c8 -> cd(2);
	f_fal_neg[no_emiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_emiss_cuts-1 ; j++) f_fal_neg[j] -> Draw("samePE");
	leg_false_neg -> Draw("same");
	el2 -> Draw("same");
	c8 -> cd(3);
	f_pos_neg[no_emiss_cuts-1] -> Draw("APE");
	for(int j = 0 ; j < no_emiss_cuts-1 ; j++) f_pos_neg[j] -> Draw("samePE");
	leg_pos_neg -> Draw("same");
	l_f_pos_neg -> Draw("same");
	el3 -> Draw("same");
	c8 -> cd(4);

	TLatex * tex_00 = new TLatex(0.05,0.95,"Begin with a sample that passes e-kinematics cuts"       );       tex_00 -> Draw(      );
	TLatex * tex_01 = new TLatex(0.05,0.85,"Denominator:"                                            );	tex_01 -> Draw("same");
	TLatex * tex_02 = new TLatex(0.05,0.80,"smeared protons pass test P_{miss}, E_{miss} cuts"       );	tex_02 -> Draw("same");
	TLatex * tex_03 = new TLatex(0.05,0.60,"False positive numerator:"                               );	tex_03 -> Draw("same");
	TLatex * tex_04 = new TLatex(0.05,0.55,"Smeared protons pass test P_{miss}, E_{miss} cuts, and"  );	tex_04 -> Draw("same");
	TLatex * tex_05 = new TLatex(0.05,0.50,"un-smeared fail P_{miss}, E_{miss} fixed"                );	tex_05 -> Draw("same");
	TLatex * tex_06 = new TLatex(0.05,0.30,"False negative numerator:"                               );	tex_06 -> Draw("same");
	TLatex * tex_07 = new TLatex(0.05,0.25,"Un-smeared protons pass fixed P_{miss}, E_{miss} and"    );       tex_07 -> Draw("same");
	TLatex * tex_08 = new TLatex(0.05,0.20,"smeared fail P_{miss}, E_{miss} cuts"   		 );	tex_08 -> Draw("same");

	c8->Modified(); c8->Update();

	TCanvas * c9 = new TCanvas("c9","",1300,800);
	c9 -> Divide(2,2);
	c9 -> cd(1);   h2_thp_y      -> Draw("COL");	l_ec_ang_reg1 -> Draw("same");
	c9 -> cd(2);   h1_y_full     -> Draw(     );	h1_y_ec       -> Draw("same");	leg_ydist_eccut -> Draw("same");
	c9 -> cd(3);   h2_thp_pmpar  -> Draw("COL");	l_ec_ang_reg2 -> Draw("same");
	c9 -> cd(4);   h1_pmpar_full -> Draw(     );	h1_pmpar_ec   -> Draw("same");	leg_ydist_eccut -> Draw("same");
	c9 -> Modified(); c9->Update();

	TCanvas * c10 = new TCanvas("c10","",1300,800);
	c10 -> Divide(2,2);
	c10 -> cd(1);	h2_y_xb       -> Draw("COL");
	c10 -> cd(2);	h1_xb_full    -> Draw(     );	h1_xb_ec -> Draw("same");	leg_ydist_eccut -> Draw("same");
	c10 -> cd(3);	h2_nu_Q2_full -> Draw("COL");	l_Q2_Nu  -> Draw("same");
	c10 -> cd(4);	h2_nu_Q2_ec   -> Draw("COL");	l_Q2_Nu  -> Draw("same");
	c10 -> Modified(); c10->Update();

	TCanvas * c11 = new TCanvas("c11","",1300,800);
	c11 -> Divide(3,2);
	c11 -> cd(1);	h2_pe_ppsm  -> Draw("COL");	l_pe_p -> Draw("same");
	c11 -> cd(2);   h2_pe_pn    -> Draw("COL");	l_pe_n -> Draw("same");
	c11 -> cd(3);	h2_pe_pp    -> Draw("COL");	l_pe_p -> Draw("same");
	c11 -> cd(4);	h1_phi_epsm -> Draw();
	c11 -> cd(5);	h1_phi_en   -> Draw();
	c11 -> cd(6);	h1_phi_ep   -> Draw();
	c11 -> Modified(); c11->Update();

	TCanvas * c12 = new TCanvas("c12","",1300,800);
	c12 -> Divide(2,2);
	c12 -> cd(1);	h1_theta_e_ps -> Draw();	h1_theta_e_n -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c12 -> cd(2);	h1_theta_ps   -> Draw();	h1_theta_n   -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c12 -> cd(3);	h2_thpq_ppq   -> Draw("COL");	l_th_pq      -> Draw("same");
	c12 -> cd(4);	h2_thnq_pnq   -> Draw("COL");	l_th_pq      -> Draw("same");
	c12 -> Modified(); c12->Update();

	TCanvas * c13 = new TCanvas("c13","",1300,800);
	c13 -> Divide(3,2);
	c13 -> cd(1);	h1_p_e_ps  -> Draw();	h1_p_e_n  -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> cd(2);	h1_th_e_ps -> Draw();	h1_th_e_n -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> cd(3);	h1_Q2_e_ps -> Draw();	h1_Q2_e_n -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> cd(4);	h1_Nu_e_ps -> Draw();	h1_Nu_e_n -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> cd(5);	h1_y_e_ps  -> Draw();	h1_y_e_n  -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> cd(6);	h1_Xb_e_ps -> Draw();	h1_Xb_e_n -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c13 -> Modified(); c13->Update();

	TCanvas * c14 = new TCanvas("c14","",1300,800);
	c14 -> Divide(2,2);
	c14 -> cd(1);	h1_p_ps   -> Draw();	h1_p_n    -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c14 -> cd(2);	h1_th_ps  -> Draw();	h1_th_n   -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c14 -> cd(3);	h1_pq_ps  -> Draw();	h1_pq_n   -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c14 -> cd(4);	h1_phi_ps -> Draw();	h1_phi_n  -> Draw("same");	leg_theta_meas_calc -> Draw("same");
	c14 -> Modified(); c14->Update();

	// --------------------------
	TCanvas * c15 = new TCanvas("c15","",1300,800);
	c15 -> Divide(2,2);
	c15 -> cd(1);

	double g_model_ratio_x[2] = {N+Z-2,N+Z+2};
	double g_model_ratio_y[2] = { 2.51, 2.51};
	double g_model_ratio_e[2] = { 0.22, 0.22};
	prep_pink_band( nu_eep_een_GK05/nu_eep , nu_eep_een_AMT/nu_eep , nu_eep_een_MMD /nu_eep , g_model_ratio_y , g_model_ratio_e );
	TGraphErrors * g_model_ratio = new TGraphErrors(2,g_model_ratio_x,g_model_ratio_y,0,g_model_ratio_e);
	g_model_ratio -> SetFillStyle(3008);
	g_model_ratio -> SetFillColor(6);
	g_model_ratio -> GetYaxis() -> SetTitle("#frac{#(e,e'p)/Z}{#(e,e'n)/N}");
	g_model_ratio -> GetXaxis() -> SetTitle("A");
	g_model_ratio -> SetTitle(tab_targ);
	g_model_ratio -> GetXaxis() -> SetRangeUser(g_model_ratio_x[0],g_model_ratio_x[1]);
	g_model_ratio -> SetMinimum(1.0);
	g_model_ratio -> SetMaximum(3.5);
	pretty_TGraphErrors(g_model_ratio);
	g_model_ratio -> Draw("AE3");

	double A_x = N+Z;
	double eep_een      = (nu_eep/Z)/(nu_een/N);
	double e_sta_eep_een = N/Z*sqrt( pow(sqrt(nu_eep)/nu_een,2) + pow(nu_eep*sqrt(nu_een)/pow(nu_een,2),2) );	// statistical error
	double e_def_eep_een = 0;//sqrt( pow(n_eff_p0*n_eff_m_e,2)+pow(n_eff_p0_e*n_eff_m,2) );				// neutron detection eff error

	// Event selection errors (Calculated with the sensitivity code and manually added here)
	double e_evt_eep_een;
	if     (tab_targ == "4He")	e_evt_eep_een = 0.;
	else if(tab_targ == "3He")	e_evt_eep_een = 0.;
	else if(tab_targ == "12C")	e_evt_eep_een = 0.;
	else	e_evt_eep_een = 0;

	double e_tot_eep_een = sqrt(e_sta_eep_een*e_sta_eep_een+e_def_eep_een*e_def_eep_een+
			e_evt_eep_een*e_evt_eep_een*eep_een*eep_een);

	cout << "**************\n[#(e,e'p)/Z]/[#(e,e'n)/N] = " << eep_een << " +/- " << e_tot_eep_een; 
	cout << "\nstatistical error: " << e_sta_eep_een;
	cout << "\nn det eff   error: " << e_def_eep_een;
	cout << "\nevt select  error: " << e_evt_eep_een*eep_een << endl;
	cout << "**************" << endl;

	cout << "Number of protons  contributing to this ratio: " << raw_nu_p << endl;
	cout << "Number of neutrons contributing to this ratio: " << raw_nu_n << endl << endl;

	TGraphErrors * g_data_ratio = new TGraphErrors(1,&A_x,&eep_een,0,&e_tot_eep_een);
	g_data_ratio -> SetMarkerStyle(20);
	g_data_ratio -> Draw("sameP");

	// ------------------------------------------------------
	// Momentum-dependent ratio:
	// Calculation
	c15 -> cd(2);

	float calc_eep_een_pdep_mom[4] = {1.50,1.73,1.97,2.21};
	float calc_eep_een_pdep_mer[4] = {0.10,0.10,0.10,0.10};
	float calc_eep_een_pdep_rat[4] = {2.55,2.54,2.50,2.45};
	float calc_eep_een_pdep_rer[4] = {0.21,0.21,0.20,0.20};
	TGraphErrors * g_calc_ratio_pdep = new TGraphErrors(4,calc_eep_een_pdep_mom,calc_eep_een_pdep_rat,calc_eep_een_pdep_mer,calc_eep_een_pdep_rer);
	g_calc_ratio_pdep -> SetFillColor(18);
	g_calc_ratio_pdep -> SetMinimum(0);
	g_calc_ratio_pdep -> SetMaximum(4);
	g_calc_ratio_pdep -> GetYaxis() -> SetTitle("#frac{#(e,e'p)/Z}{#(e,e'n)/N}");
	g_calc_ratio_pdep -> GetXaxis() -> SetTitle("p_{N} [GeV]");
	g_calc_ratio_pdep -> SetTitle(tab_targ);
	pretty_TGraphErrors(g_calc_ratio_pdep);
	g_calc_ratio_pdep -> Draw("a2");
	// Data
	double eep_een_pdep      [bins_eepn];
	double e_sta_eep_een_pdep[bins_eepn];
	double e_tot_eep_een_pdep[bins_eepn];
	for(int i1 = 0 ; i1 < bins_eepn   ; i1++) {
		eep_een_pdep      [i1] = (nu_eep_p[i1]/Z)/(nu_een_p[i1]/N);
		e_sta_eep_een_pdep[i1] = N/Z*sqrt( pow(sqrt(nu_eep_p[i1])/nu_een_p[i1],2) + pow(nu_eep_p[i1]*sqrt(nu_een_p[i1])/pow(nu_een_p[i1],2),2) );
		e_tot_eep_een_pdep[i1] = sqrt(e_sta_eep_een_pdep[i1]*e_sta_eep_een_pdep[i1]+e_def_eep_een*e_def_eep_een+
				e_evt_eep_een*e_evt_eep_een*eep_een*eep_een);
	}
	TGraphErrors * g_data_ratio_pdep = new TGraphErrors(bins_eepn,p_ctr_eepn,eep_een_pdep,p_err_eepn,e_tot_eep_een_pdep);
	g_data_ratio_pdep -> Draw("sameP");

	c15 -> Modified(); c15->Update();

	// ------------------------------------------------------
	// Additional plots
	TCanvas * c16 = new TCanvas("c16","",1300,800);
	c16 -> Divide(2,2);
	c16 -> cd(1);   h2_p_pm_em   -> Draw("COLZ");
	c16 -> cd(2);   h2_n_pm_em   -> Draw("COLZ");
	c16 -> cd(3);   h2_psm_pm_em -> Draw("COLZ");
	c16 -> Modified(); c16->Update();

	TCanvas * c17 = new TCanvas("c17","",1300,800);
	c17 -> Divide(2,2);
	c17 -> cd(1);	h2_p_pm_em_1   -> Draw("COLZ");
	c17 -> cd(2);	h2_n_pm_em_1   -> Draw("COLZ");
	c17 -> cd(3);	h2_psm_pm_em_1 -> Draw("COLZ");
	c17 -> Modified(); c17->Update();

	TCanvas * c18 = new TCanvas("c18","",1300,800);
	c18 -> Divide(2,2);
	c18 -> cd(1);	h2_pm_pmsm -> Draw("COLZ");
	c18 -> cd(2);	h2_em_emsm -> Draw("COLZ");
	c18 -> cd(3);	h2_ph_phsm -> Draw("COLZ");
	c18 -> cd(4);	h2_ph_pmsm -> Draw("COLZ");
	c18 -> Modified(); c18->Update();

	TCanvas * c19 = new TCanvas("c19","",1300,800);
	c19 -> Divide(2,2);
	c19 -> cd(1);	h2_pm_pmsm_pm_diff -> Draw("COLZ");
	c19 -> cd(2);	h2_em_emsm_em_diff -> Draw("COLZ");
	c19 -> cd(3);	h2_ph_phsm_ph_diff -> Draw("COLZ");
	c19 -> cd(4);	h2_ph_pmsm_ph_diff -> Draw("COLZ");
	c19 -> Modified(); c19->Update();

	TCanvas * c20 = new TCanvas("c20","",1300,800);
        c20 -> Divide(2,2);
        c20 -> cd(1);	gPad -> SetLogz();	h2_but_pm_em_0 -> Draw("COLZ");
	c20 -> cd(2);	gPad -> SetLogz();	h2_but_pm_em_1 -> Draw("COLZ");
	c20 -> cd(3);	gPad -> SetLogz();	h2_but_pm_ph   -> Draw("COLZ");
	c20 -> cd(4);	gPad -> SetLogz();	h2_but_em_ph   -> Draw("COLZ");
	c20 -> Modified(); c20->Update();

	TCanvas * c21 = new TCanvas("c21","",1300,800);
	c21 -> Divide(2,2);
        c21 -> cd(1);    h2_p_Nu_y_QEcut    -> Draw("COL");      b_Nu_y -> Draw("same");
        c21 -> cd(2);    h2_p_Nu_y_QEcut2   -> Draw("COL");      b_Nu_y -> Draw("same");
	c21 -> cd(3);    h2_p_Q2_thpq_QEcut -> Draw("COL");      b_Q2_th_pq -> Draw("same");
	c21 -> cd(4);    h2_p_Q2_thpq_QEcut2-> Draw("COL");      b_Q2_th_pq -> Draw("same");
	c21 -> Modified(); c21->Update();

	// --------------------------
	// Printing plots to pdf
	c4  ->Print("./Results/plots_MF_"+tab_targ+".pdf(");
	c5  ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c7  ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c8  ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c9  ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c10 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c11 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c12 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c13 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c14 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );	
	c15 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c16 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c17 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c18 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c19 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c20 ->Print("./Results/plots_MF_"+tab_targ+".pdf" );
	c21 ->Print("./Results/plots_MF_"+tab_targ+".pdf)");

	// --------------------------
	// Creating output table
	ofstream output_table;
	output_table.open("./Results/table_MF_"+tab_targ+".dat");
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
double y_scale_func(TLorentzVector y_q, double y_Mtar, double y_Mres, double y_Mnuc){
	// See equation 9a from "Scaling in inclusive electron-nucleus scattering
	// D.B. Day, J.S. McCarthy, T.W. Donnelly, and I. Sick
	// Annu. Rev. Nucl. Part. Sci 1990. 40: 357-409
	double q3       = y_q.Rho();
	double nu       = y_q.E  ();
	double W        = sqrt( pow(y_Mtar+nu,2) - pow(q3,2) );
	double Lambda   = ( pow(y_Mres,2) - pow(y_Mnuc,2) + pow(W,2) )/2.;
	return ( (y_Mtar+nu)*sqrt( pow(Lambda,2) - pow(y_Mres*W,2) ) - q3*Lambda )/pow(W,2);     
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

