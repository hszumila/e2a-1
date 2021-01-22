#include "TH1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TROOT.h"
#include <vector>
#include "unistd.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TCut.h"
#include "Riostream.h"
#include "TNtupleD.h"
#include "TCut.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPaveText.h"
#include <iostream>

#include "Fiducial.h"

using namespace std;


const double Me      = .000511;         // Electron     mass in GeV
const double Mp      = 0.93827;         // Proton       mass in GeV
const double Mn      = 0.93957;         // Neutron      mass in GeV
const double Mpip    = 0.13957;         // Charged pion mass in GeV
const double He3_be  = 0.00770;         // 3He binding energy in GeV
const double M_tar   = 2*Mp+Mn-He3_be;  // 3Helium  mass in GeV (target)
const double c_cm_s  = 2.99792E+10; //speed of the light in cm/s
const double ns2s    = 1.0E-09 ;
const double corr_pp[3] = {0.6621,0.7795,0.2584};//bg subtraction for e,e'pp
const double corr_pppipi[2] = {0.7008,0.5480};//bg subtraction for e,e'pppipi

double max_outta_three ( double A , double B, double C );
Double_t func_fit(Double_t *x, Double_t *par);
void fit_gfunction(TH1F * gPlot, float minx, float maxx);
void fit_gfunctionExt(TH1F * gPlot, float minx, float maxx);

Double_t gaussianPeak(Double_t *x, Double_t *par);
Double_t background(Double_t *x, Double_t *par);
Double_t fitSB(Double_t *x, Double_t *par);

bool p_points_to_EC_fid(TVector3 pm, double vtx_z, double dist){
	double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;
	double pi = 22./7.;

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


int main(){

  //  option = 0: 4.4 GeV data
  //         = 1: 2.2 GeV data
  int option = 1;
  
  bool eppn, eppnpipi;
  
  //Define Branch variables and other variables needed to determine parameters
  int nEntries;
  
  const double Mn_up   = 1.0000;		// Delta neutron mass in GeV for mass distribution cut
  const double Mn_do   = 0.8000;		// Delta neutron mass in GeV for mass distribution cut
  
  double Ebeam, M_miss_eppn, M_miss_eppnpipi, nBeta, Rcalc, cos_th_miss_ec_eppn, cos_th_miss_ec_eppnpipi;
  TVector3 EC_hit_pos, P_miss_eppn, P_miss_eppnpipi, u1;
  TLorentzVector mom_e,mom_p1,mom_p2, E0, V4_tar, mom_pip, mom_pim;
  double sig_e_p1, sig_e_p2, sig_p1p2, sig_e_pip, sig_e_pim, sig_pip_pim, temp1, temp2,
    Z_e_p1_abs, Z_e_p2_abs, Z_p1p2_abs, Z_e_pip_abs, Z_e_pim_abs, Z_pip_pim_abs, vz_o_sig;
  // --------------------------------------------------------------------------
  //Read in data files
  TFile *f;
  
  if(option==0){
    f = new TFile("../../../data/n_eff/skim_He3_pp.root");
    Ebeam = 4.4;
  }
  else if(option==1){
    f = new TFile("../../../data/n_eff/skim_He3_pp_2261.root");
    Ebeam = 2.2;
  }
  TTree * t = (TTree*)f->Get("T");
  TString outputpdf = Form("out_neff/output_%.1f",Ebeam);
  // --------------------------------------------------------------------------
  int nParticles, nProtons, nNeutrons, nPiplus, nPiminus;
  int Part_type [20], stat_dc  [20], stat_sc[20], stat_ec[20];
  double t0, QSq, xB, Nu;
  double mom_x  [20], mom_y    [20], mom_z  [20], Mass   [20], charge[20], beta[20], vtx_z_cor[20],
    ec_time[20], ec_path  [20], ec_in  [20], ec_out [20], ec_tot[20], ec_x[20], ec_y     [20], ec_z[20],
    ec_u   [20], ec_v     [20], ec_w   [20];
  //Define Branch Addresses
  t->SetBranchAddress("nParticles",&nParticles);
  t->SetBranchAddress("nProtons"  ,&nProtons  );
  t->SetBranchAddress("nNeutrons" ,&nNeutrons );
  t->SetBranchAddress("nPiplus" ,&nPiplus );
  t->SetBranchAddress("nPiminus" ,&nPiminus );
  t->SetBranchAddress("Q2",&QSq);
  t->SetBranchAddress("Xb",&xB);
  t->SetBranchAddress("Nu",&Nu);
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

  TH1F * h1_pmiss_epp_corr        = new TH1F("h1_pmiss_epp_corr"        , "^{3}He(e,e'pp)n;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
  TH1F * h1_pmiss_epppipi_corr    = new TH1F("h1_pmiss_epppipi_corr"    , "^{3}He(e,e'pp#pi^{+}#pi^{-})n;p_{miss} [GeV];Counts"   , 20 , 0.5 , 2.1 );
  TH1F * h1_pmiss_both_inc_corr   = new TH1F("h1_pmiss_both_inc_corr"   , "both, inclusive;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
  
  TH1F * h1_pmiss_eppn       = new TH1F("h1_pmiss_eppn"       , "^{3}He(e,e'ppn);p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
  TH1F * h1_pmiss_epppipin   = new TH1F("h1_pmiss_epppipin"   , "^{3}He(e,e'pp#pi^{+}#pi^{-}n);p_{miss} [GeV];Counts"   , 20 , 0.5 , 2.1 );
  TH1F * h1_pmiss_both_exc   = new TH1F("h1_pmiss_both_exc"   , "both, exclusive;p_{miss} [GeV];Counts"                 , 20 , 0.5 , 2.1 );
  
  TH1F * h1_det_eff_eppn     = new TH1F("h1_det_eff_eppn"     , ""                                                      , 20 , 0.5 , 2.1 );
  TH1F * h1_det_eff_eppnpipi = new TH1F("h1_det_eff_eppnpipi" , ""                                                      , 20 , 0.5 , 2.1 );
  TH1F * h1_det_eff_both     = new TH1F("h1_det_eff_both"     , ""                                                      , 20 , 0.5 , 2.1 );

  TH1F * h1_det_eff_eppn_corr     = new TH1F("h1_det_eff_eppn_corr"     , ""                                                      , 20 , 0.5 , 2.1 );
  TH1F * h1_det_eff_eppnpipi_corr = new TH1F("h1_det_eff_eppnpipi_corr" , ""                                                      , 20 , 0.5 , 2.1 );
  TH1F * h1_det_eff_both_corr     = new TH1F("h1_det_eff_both_corr"     , ""                                                      , 20 , 0.5 , 2.1 );


  TH1F *h1_Mmiss_eppn_inc0 = new TH1F("h1_Mmiss_eppn_inc0","0.5<p_{miss}<1 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppn_inc1 = new TH1F("h1_Mmiss_eppn_inc1","1<p_{miss}<1.5 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppn_inc2 = new TH1F("h1_Mmiss_eppn_inc2","1.5<p_{miss}<2.1 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppn_exc0 = new TH1F("h1_Mmiss_eppn_exc0","0.5<p_{miss}<1 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppn_exc1 = new TH1F("h1_Mmiss_eppn_exc1","1<p_{miss}<1.5 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppn_exc2 = new TH1F("h1_Mmiss_eppn_exc2","1.5<p_{miss}<2.1 GeV;m_{miss} [GeV]",30,0.4,1.4);

  TH1F *h1_Mmiss_eppin_inc0 = new TH1F("h1_Mmiss_eppin_inc0","0.5<p_{miss}<1.2 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppin_inc1 = new TH1F("h1_Mmiss_eppin_inc1","1.2<p_{miss}<2.1 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppin_exc0 = new TH1F("h1_Mmiss_eppin_exc0","0.5<p_{miss}<1.2 GeV;m_{miss} [GeV]",30,0.4,1.4);
  TH1F *h1_Mmiss_eppin_exc1 = new TH1F("h1_Mmiss_eppin_exc1","1.2<p_{miss}<2.1 GeV;m_{miss} [GeV]",30,0.4,1.4);

  //TH2D * h1_xy_eppn= new TH2D("h1_xy_eppn","^{3}He(e,e'pp)n;EC_{x} [cm];EC_{y} [cm];Counts"                 ,100,-500.,500.,100,-500,500);
  //TH2D * h1_xy_eppnpipi= new TH2D("h1_xy_eppnpipi","^{3}He(e,e'pp)n;EC_{x} [cm];EC_{y} [cm];Counts"                 ,100,-500.,500.,100,-500,500);
  
  cout << "Histograms have been defined" << endl;

  Fiducial fid_params(4461,2250,5996,"3He",true);

  
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
  c0->Print(outputpdf+".pdf(");

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

    if (!fid_params.pFiducialCut(P_miss_eppn)) {continue;}
    //if(P_miss_eppn.Theta()>44.*3.14/180.) continue;

    // --------------------------------------------------------------------------
    // Cuts
    // -----------------------------------
    // 3He(e,e'ppn)
    // -----------------------------------
    //fill exclusive bg plot:
    Z_e_p1_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p1])/sig_e_p1;
    Z_e_p2_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p2])/sig_e_p2;
    Z_p1p2_abs = TMath::Abs(vtx_z_cor[p1]-vtx_z_cor[p2])/sig_p1p2;
    vz_o_sig = max_outta_three(Z_e_p1_abs,Z_e_p2_abs,Z_p1p2_abs);
    if(!p_points_to_EC_fid(P_miss_eppn,vtx_z_cor[ 0],0.)) {continue;}

    //hole in acceptance
    /*
    double phi=P_miss_eppn.Phi();
    double theta = P_miss_eppn.Theta()*180./3.14;  
    if (phi < -M_PI/6.) phi+= 2.*M_PI;
    phi *= 180./M_PI; 
    if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
    if (phi>150 && phi<200 && theta>35 && theta<40){continue;}
    */
    
    if(vz_o_sig<3.){
      //fill inclusive mmiss plots, with vertex cut
      if (P_miss_eppn.Mag()>0.5 && P_miss_eppn.Mag()<=1.0){h1_Mmiss_eppn_inc0->Fill(M_miss_eppn);}
      if (P_miss_eppn.Mag()>1.0 && P_miss_eppn.Mag()<=1.5){h1_Mmiss_eppn_inc1->Fill(M_miss_eppn);}
      if (P_miss_eppn.Mag()>1.5 && P_miss_eppn.Mag()<=2.1){h1_Mmiss_eppn_inc2->Fill(M_miss_eppn);}
    }
    
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

	  //corrected the inclusive efficiency
	  if (P_miss_eppn.Mag()<=1.){
	    h1_pmiss_epp_corr        -> Fill(P_miss_eppn    .Mag(),corr_pp[0]);
	    h1_pmiss_both_inc_corr   -> Fill(P_miss_eppn    .Mag(),corr_pp[0]);
	  }
	  if (P_miss_eppn.Mag()>1 && P_miss_eppn.Mag()<=1.5){
	    h1_pmiss_epp_corr        -> Fill(P_miss_eppn    .Mag(),corr_pp[1]);
	    h1_pmiss_both_inc_corr   -> Fill(P_miss_eppn    .Mag(),corr_pp[1]);
	  }
	  if (P_miss_eppn.Mag()>1.5 && P_miss_eppn.Mag()<=2.1){
	    h1_pmiss_epp_corr        -> Fill(P_miss_eppn    .Mag(),corr_pp[2]);
	    h1_pmiss_both_inc_corr   -> Fill(P_miss_eppn    .Mag(),corr_pp[2]);
	  }

	  
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
		      h1_det_eff_eppn_corr   -> Fill(P_miss_eppn.Mag());
		      h1_det_eff_both_corr   -> Fill(P_miss_eppn.Mag());
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    //calc Mmiss spectra with no Mmiss cut for bg analysis:
    cos_th_miss_ec_eppn = P_miss_eppn*EC_hit_pos/P_miss_eppn.Mag()/EC_hit_pos.Mag();
    if(vz_o_sig<3. && P_miss_eppn.Mag()> 0.5 && n1!=0 && stat_ec[n1]>0 && stat_dc[n1]<=0 && stat_sc[n1]<=0 &&beta[n1]<0.95 && ec_tot[n1]>0 && ec_time[n1]>0 && cos_th_miss_ec_eppn>0.995 && ec_u[n1]>40 && ec_v[n1]<360 && ec_w[n1]<395) {
	if (P_miss_eppn.Mag()>0.5 && P_miss_eppn.Mag()<=1.0){h1_Mmiss_eppn_exc0->Fill(M_miss_eppn);}
	if (P_miss_eppn.Mag()>1.0 && P_miss_eppn.Mag()<=1.5){h1_Mmiss_eppn_exc1->Fill(M_miss_eppn);}
	if (P_miss_eppn.Mag()>1.5 && P_miss_eppn.Mag()<=2.1){h1_Mmiss_eppn_exc2->Fill(M_miss_eppn);}
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
    if (!fid_params.pFiducialCut(P_miss_eppnpipi)) {continue;}

    // --------------------------------------------------------------------------
    // Cuts
    // -----------------------------------
    //fill exclusive bg plot:
    Z_e_p1_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p1])/sig_e_p1   ;
    Z_e_p2_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p2])/sig_e_p2   ;
    Z_p1p2_abs    = TMath::Abs(vtx_z_cor[ p1]-vtx_z_cor[ p2])/sig_p1p2   ;
    Z_e_pip_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pip])/sig_e_pip  ;
    Z_e_pim_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pim])/sig_e_pim  ;
    Z_pip_pim_abs = TMath::Abs(vtx_z_cor[pip]-vtx_z_cor[pim])/sig_pip_pim;
    temp1    = max_outta_three(Z_e_p1_abs ,Z_e_p2_abs ,Z_p1p2_abs   );
    temp2    = max_outta_three(Z_e_pip_abs,Z_e_pim_abs,Z_pip_pim_abs);
    vz_o_sig = TMath::Max(temp1,temp2);
    if(!p_points_to_EC_fid(P_miss_eppnpipi,vtx_z_cor[ 0],0.)) {continue;}

    if(vz_o_sig<3.){
      //fill inclusive mmiss plots, with vertex cut
      if (P_miss_eppnpipi.Mag()>0.5 && P_miss_eppnpipi.Mag()<=1.2){h1_Mmiss_eppin_inc0->Fill(M_miss_eppnpipi);}
      if (P_miss_eppnpipi.Mag()>1.2 && P_miss_eppnpipi.Mag()<=2.1){h1_Mmiss_eppin_inc1->Fill(M_miss_eppnpipi);}
    }
    
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


	  //corrected the inclusive efficiency
	  if (P_miss_eppnpipi.Mag()<=1.2){
	    h1_pmiss_epppipi_corr        -> Fill(P_miss_eppnpipi    .Mag(),corr_pppipi[0]);
	    h1_pmiss_both_inc_corr   -> Fill(P_miss_eppnpipi    .Mag(),corr_pppipi[0]);
	  }
	  if (P_miss_eppnpipi.Mag()>1.2 && P_miss_eppnpipi.Mag()<=2.1){
	    h1_pmiss_epppipi_corr        -> Fill(P_miss_eppnpipi    .Mag(),corr_pppipi[1]);
	    h1_pmiss_both_inc_corr   -> Fill(P_miss_eppnpipi    .Mag(),corr_pppipi[1]);
	  }
	  	  
	  
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
		      h1_det_eff_eppnpipi_corr -> Fill(P_miss_eppnpipi.Mag());
		      h1_det_eff_both_corr     -> Fill(P_miss_eppnpipi.Mag());
		    }
		  }
		}
	      }	
	    }
	  }
	}
      }
    }
    //calc Mmiss spectra with no Mmiss cut for bg analysis:
    cos_th_miss_ec_eppnpipi = P_miss_eppnpipi*EC_hit_pos/P_miss_eppnpipi.Mag()/EC_hit_pos.Mag();
    if(vz_o_sig<3. && n1!=0 && stat_ec[n1]>0 && stat_dc[n1]<=0 && stat_sc[n1]<=0 &&beta[n1]<0.95 && ec_tot[n1]>0 && ec_time[n1]>0 && cos_th_miss_ec_eppnpipi>0.995 && ec_u[n1]>40 && ec_v[n1]<360 && ec_w[n1]<395) {
      if (P_miss_eppnpipi.Mag()>0.5 && P_miss_eppnpipi.Mag()<=1.2){h1_Mmiss_eppin_exc0->Fill(M_miss_eppnpipi);}
      if (P_miss_eppnpipi.Mag()>1.2 && P_miss_eppnpipi.Mag()<=2.1){h1_Mmiss_eppin_exc1->Fill(M_miss_eppnpipi);}
    }
    // -------------------------------------------------------------------------- 
  }
  cout << "Done with the loop" << endl;
  
  // --------------------------------------------------------------------------
  //Calculate the signal to background
  //e,e'pp reaction
  TF1 *fitFcn_pin0 = new TF1("fitFcn_pin0",fitSB,0.4,1.3,7);
  fitFcn_pin0->SetParameters(100,-300,300,-300,150,0.9,0.05);
  fitFcn_pin0->SetParLimits(5,0.8,0.94);
  h1_Mmiss_eppn_inc0->Fit("fitFcn_pin0","R");
  TF1 *fSignal0 = new TF1("fSignal0",gaussianPeak,0.4,1.3,3);
  fSignal0->SetLineColor(kBlue);
  TF1 *fBackground0 = new TF1("fBackground0",background,0.4,1.3,4);
  Double_t param0[7];
  fitFcn_pin0->GetParameters(param0);
  fSignal0->SetParameters(&param0[4]);
  fBackground0->SetParameters(&param0[0]);

  TF1 *fitFcn_pin1 = new TF1("fitFcn_pin1",fitSB,0.8,1.3,7);
  fitFcn_pin1->SetParameters(100,-300,300,-300,150,0.9,0.05);
  fitFcn_pin1->SetParLimits(5,0.85,0.94);
  h1_Mmiss_eppn_inc1->Fit("fitFcn_pin1","R");
  TF1 *fSignal1 = new TF1("fSignal1",gaussianPeak,0.8,1,3);
  fSignal1->SetLineColor(kBlue);
  TF1 *fBackground1 = new TF1("fBackground1",background,0.8,1.3,4);
  Double_t param1[7];
  fitFcn_pin1->GetParameters(param1);
  fSignal1->SetParameters(&param1[4]);
  fBackground1->SetParameters(&param1[0]);
  
  TF1 *fitFcn_pin2 = new TF1("fitFcn_pin2",fitSB,0.4,1.3,7);
  fitFcn_pin2->SetParameters(100,-300,300,-300,150,0.9,0.05);
  fitFcn_pin2->SetParLimits(5,0.8,0.95);
  h1_Mmiss_eppn_inc2->Fit("fitFcn_pin2","R");
  TF1 *fSignal2 = new TF1("fSignal2",gaussianPeak,0.4,1.3,3);
  fSignal2->SetLineColor(kBlue);
  TF1 *fBackground2 = new TF1("fBackground2",background,0.4,1.3,4);
  Double_t param2[7];
  fitFcn_pin2->GetParameters(param2);
  fSignal2->SetParameters(&param2[4]);
  fBackground2->SetParameters(&param2[0]);
  
  double sbsum[3],sig[3];

  sbsum[0]=fitFcn_pin0->Integral(Mn_do,Mn_up);
  sbsum[1]=fitFcn_pin1->Integral(Mn_do,Mn_up);
  sbsum[2]=fitFcn_pin2->Integral(Mn_do,Mn_up);

  sig[0]=fSignal0->Integral(Mn_do,Mn_up);
  sig[1]=fSignal1->Integral(Mn_do,Mn_up);
  sig[2]=fSignal2->Integral(Mn_do,Mn_up);

  cout<<"S: "<<sig[0]<<" S+B: "<<sbsum[0]<<" ratio: "<<sig[0]/sbsum[0]<<endl;
  cout<<"S: "<<sig[1]<<" S+B: "<<sbsum[1]<<" ratio: "<<sig[1]/sbsum[1]<<endl;
  cout<<"S: "<<sig[2]<<" S+B: "<<sbsum[2]<<" ratio: "<<sig[2]/sbsum[2]<<endl;

  
  TCanvas * c4 = new TCanvas("c4","3He(e,e'pp)n",1200,900);
  c4 -> Divide(3,2); 
  c4 -> cd(1);	h1_Mmiss_eppn_inc0     -> Draw();
  fSignal0->SetLineColor(kBlue); fBackground0->SetLineColor(kGreen);
  fSignal0->Draw("same"); fBackground0->Draw("same");
  TPaveText *pt0 = new TPaveText(.5,50,0.75,100);
  pt0->AddText(Form("S/(S+B) = %.4f",sig[0]/sbsum[0]));
  pt0->Draw("same");
  c4 -> cd(2);	h1_Mmiss_eppn_inc1     -> Draw();
  fSignal1->SetLineColor(kBlue); fBackground1->SetLineColor(kGreen);
  fSignal1->Draw("same"); fBackground1->Draw("same");
  TPaveText *pt1 = new TPaveText(.5,15,0.75,25);
  pt1->AddText(Form("S/(S+B) = %.4f",sig[1]/sbsum[1]));
  pt1->Draw("same");
  c4 -> cd(3);	h1_Mmiss_eppn_inc2     -> Draw();
  fSignal2->SetLineColor(kBlue); fBackground2->SetLineColor(kGreen);
  fSignal2->Draw("same"); fBackground2->Draw("same");
  TPaveText *pt2 = new TPaveText(.5,10,0.75,20);
  pt2->AddText(Form("S/(S+B) = %.4f",sig[2]/sbsum[2]));
  pt2->Draw("same");
  c4 -> cd(4);	h1_Mmiss_eppn_exc0     -> Draw();
  c4 -> cd(5);	h1_Mmiss_eppn_exc1     -> Draw();
  c4 -> cd(6);	h1_Mmiss_eppn_exc2     -> Draw();
  c4->Print(outputpdf+".pdf(");
  
  //e,e'pppippim reaction _ppin
  TF1 *fitFcn_ppin0 = new TF1("fitFcn_ppin0",fitSB,0.7,1.3,7);
  fitFcn_ppin0->SetParameters(100,-300,300,-300,150,0.9,0.05);
  fitFcn_ppin0->SetParLimits(5,0.88,0.95);
  h1_Mmiss_eppin_inc0->Fit("fitFcn_ppin0","R");
  TF1 *fSignal0pi = new TF1("fSignal0pi",gaussianPeak,0.7,1.3,3);
  fSignal0pi->SetLineColor(kBlue);
  TF1 *fBackground0pi = new TF1("fBackground0pi",background,0.7,1.3,4);
  Double_t param0pi[7];
  fitFcn_ppin0->GetParameters(param0pi);
  fSignal0pi->SetParameters(&param0pi[4]);
  fBackground0pi->SetParameters(&param0pi[0]);
  
  TF1 *fitFcn_ppin1 = new TF1("fitFcn_ppin1",fitSB,0.4,1.3,7);
  fitFcn_ppin1->SetParameters(100,-300,300,-300,150,0.9,0.05);
  fitFcn_ppin1->SetParLimits(5,0.88,0.95);
  h1_Mmiss_eppin_inc1->Fit("fitFcn_ppin1","R");
  TF1 *fSignal1pi = new TF1("fSignal1pi",gaussianPeak,0.88,0.97,3);
  fSignal1pi->SetLineColor(kBlue);
  TF1 *fBackground1pi = new TF1("fBackground1pi",background,0.4,1.3,4);
  Double_t param1pi[7];
  fitFcn_ppin1->GetParameters(param1pi);
  fSignal1pi->SetParameters(&param1pi[4]);
  fBackground1pi->SetParameters(&param1pi[0]);
  
  
  double sbsumpi[2],sigpi[2];

  sbsumpi[0]=fitFcn_ppin0->Integral(Mn_do,Mn_up);
  sbsumpi[1]=fitFcn_ppin1->Integral(Mn_do,Mn_up);

  sigpi[0]=fSignal0pi->Integral(Mn_do,Mn_up);
  sigpi[1]=fSignal1pi->Integral(Mn_do,Mn_up);

  cout<<"S: "<<sigpi[0]<<" S+B: "<<sbsumpi[0]<<" ratio: "<<sigpi[0]/sbsumpi[0]<<endl;
  cout<<"S: "<<sigpi[1]<<" S+B: "<<sbsumpi[1]<<" ratio: "<<sigpi[1]/sbsumpi[1]<<endl;

  
  TCanvas * c5 = new TCanvas("c5","3He(e,e'pppipi)n",1200,900);
  c5 -> Divide(2,2); 
  c5 -> cd(1);	h1_Mmiss_eppin_inc0     -> Draw();
  fSignal0pi->SetLineColor(kBlue); fBackground0pi->SetLineColor(kGreen);
  fSignal0pi->Draw("same"); fBackground0pi->Draw("same");
  TPaveText *pt0pi = new TPaveText(.5,10,0.75,20);
  pt0pi->AddText(Form("S/(S+B) = %.4f",sigpi[0]/sbsumpi[0]));
  pt0pi->Draw("same");
  
  c5 -> cd(2);	h1_Mmiss_eppin_inc1     -> Draw();
  fSignal1pi->SetLineColor(kBlue); fBackground1pi->SetLineColor(kGreen);
  fSignal1pi->Draw("same"); fBackground1pi->Draw("same");
  TPaveText *pt1pi = new TPaveText(.5,5,0.75,8);
  pt1pi->AddText(Form("S/(S+B) = %.4f",sigpi[1]/sbsumpi[1]));
  pt1pi->Draw("same");
 
  c5 -> cd(3);	h1_Mmiss_eppin_exc0     -> Draw();
  c5 -> cd(4);	h1_Mmiss_eppin_exc1     -> Draw();
  c5->Print(outputpdf+".pdf(");

 
  // Calculating neutron detection efficiency
  h1_det_eff_eppn     -> Divide( h1_pmiss_epp      );
  h1_det_eff_eppnpipi -> Divide( h1_pmiss_epppipi  );
  h1_det_eff_both     -> Divide( h1_pmiss_both_inc );

  h1_det_eff_eppn_corr     -> Divide( h1_pmiss_epp_corr      );
  h1_det_eff_eppnpipi_corr -> Divide( h1_pmiss_epppipi_corr  );
  h1_det_eff_both_corr     -> Divide( h1_pmiss_both_inc_corr );
  
  double error_eppn, error_eppnpipi, error_both, A, B, eA, eB;
  
  for( int i = 1 ; i <= h1_det_eff_eppn -> GetXaxis()->GetNbins(); i++){
    A  = h1_pmiss_eppn -> GetBinContent(i);
    B  = h1_pmiss_epp  -> GetBinContent(i);
    eA = h1_pmiss_eppn -> GetBinError  (i);
    eB = h1_pmiss_epp  -> GetBinError  (i);
    error_eppn = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
    h1_det_eff_eppn -> SetBinError(i,error_eppn);

    B  = h1_pmiss_epp_corr  -> GetBinContent(i);
    eB = h1_pmiss_epp_corr  -> GetBinError  (i);
    error_eppn = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
    h1_det_eff_eppn_corr -> SetBinError(i,error_eppn);
  }
  
  for( int i = 1 ; i <= h1_det_eff_eppnpipi -> GetXaxis()->GetNbins(); i++){
    A  = h1_pmiss_epppipin -> GetBinContent(i);
    B  = h1_pmiss_epppipi  -> GetBinContent(i);
    eA = h1_pmiss_epppipin -> GetBinError  (i);
    eB = h1_pmiss_epppipi  -> GetBinError  (i);
    error_eppnpipi = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
    h1_det_eff_eppnpipi -> SetBinError(i,error_eppnpipi);
    
    B  = h1_pmiss_epppipi_corr  -> GetBinContent(i);
    eB = h1_pmiss_epppipi_corr  -> GetBinError  (i);
    error_eppnpipi = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
    h1_det_eff_eppnpipi_corr -> SetBinError(i,error_eppnpipi);
  }
  
  for( int i = 1 ; i <= h1_det_eff_both -> GetXaxis()->GetNbins(); i++){
    A  = h1_pmiss_both_exc -> GetBinContent(i);
    B  = h1_pmiss_both_inc -> GetBinContent(i);
    eA = h1_pmiss_both_exc -> GetBinError  (i);
    eB = h1_pmiss_both_inc -> GetBinError  (i);
    error_both = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
    h1_det_eff_both -> SetBinError(i,error_both);

    B  = h1_pmiss_both_inc_corr -> GetBinContent(i);
    eB = h1_pmiss_both_inc_corr -> GetBinError  (i);
    error_both = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );    
    h1_det_eff_both_corr -> SetBinError(i,error_both);
  }
  
  h1_det_eff_eppn     -> SetTitle("^{3}He(e,e'ppn)/^{3}He(e,e'pp)n"                             );
  h1_det_eff_eppn_corr     -> SetTitle("^{3}He(e,e'ppn)/^{3}He(e,e'pp)n"                             );
  h1_det_eff_eppnpipi -> SetTitle("^{3}He(e,e'pp#pi^{+}#pi^{-}n)/^{3}He(e,e'pp#pi^{+}#pi^{-})n" );
  h1_det_eff_eppnpipi_corr -> SetTitle("^{3}He(e,e'pp#pi^{+}#pi^{-}n)/^{3}He(e,e'pp#pi^{+}#pi^{-})n" );
  h1_det_eff_both     -> SetTitle("Combined neutron detection efficiency"                       );
  h1_det_eff_both_corr     -> SetTitle("Combined neutron detection efficiency, corrected"                       );
  
  h1_det_eff_eppn     -> GetYaxis() -> SetTitle("neut. det. eff.");
  h1_det_eff_eppn_corr     -> GetYaxis() -> SetTitle("neut. det. eff.");
  h1_det_eff_eppnpipi -> GetYaxis() -> SetTitle("neut. det. eff.");
  h1_det_eff_eppnpipi_corr -> GetYaxis() -> SetTitle("neut. det. eff.");
  h1_det_eff_both     -> GetYaxis() -> SetTitle("neut. det. eff.");
  h1_det_eff_both_corr     -> GetYaxis() -> SetTitle("neut. det. eff.");
  
  h1_det_eff_eppn     -> SetFillColor(10);
  h1_det_eff_eppn_corr     -> SetFillColor(10);
  h1_det_eff_eppnpipi -> SetFillColor(10);
  h1_det_eff_eppnpipi_corr -> SetFillColor(10);
  h1_det_eff_both     -> SetFillColor(10);
  h1_det_eff_both_corr     -> SetFillColor(10);
  
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
  c1->Print(outputpdf+".pdf(");

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
  c2->Print(outputpdf+".pdf(");
 
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
  c3->Print(outputpdf+".pdf(");

  TCanvas *c6 = new TCanvas("c6","Total corrected neutron efficiency",1200,900);
  c6 -> Divide(3,3);
  c6 -> cd(1);	h1_pmiss_epp_corr        -> Draw();
  c6 -> cd(2);    h1_pmiss_eppn       -> Draw();
  c6 -> cd(3);    h1_det_eff_eppn_corr     -> Draw("BARE");
  fit_gfunction( h1_det_eff_eppn_corr,0.5,2.1);
  
  c6 -> cd(4); 	h1_pmiss_epppipi_corr    -> Draw();
  c6 -> cd(5);    h1_pmiss_epppipin   -> Draw();
  c6 -> cd(6);    h1_det_eff_eppnpipi_corr -> Draw("BARE");
  fit_gfunction( h1_det_eff_eppnpipi_corr,0.5,2.1);
  
  c6 -> cd(7); 	h1_pmiss_both_inc_corr   -> Draw();
  c6 -> cd(8);	h1_pmiss_both_exc   -> Draw();
  c6 -> cd(9);	h1_det_eff_both_corr     -> Draw("BARE");
  fit_gfunction( h1_det_eff_both_corr,0.5,2.1);
  c6->Print(outputpdf+".pdf)");
  
  
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

Double_t func_fitExt(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t value;
  if     (xx <  par[0]){value = par[1]*xx    +par[2];}
  else if(xx >= par[0] && xx<par[3]){value = par[4]*xx+(par[1]-par[4])*par[0]+par[2];}
  else {value = par[4]*par[3]+(par[1]-par[4])*par[0]+par[2];}
  return value;
}
// ====================================================================================================================================================
void fit_gfunction(TH1F * gPlot, float minx, float maxx){
  Double_t par[3] = {1.5,0.3,0.1};
  TF1* Func = new TF1("Func", func_fit, 0.6,2.1,3);
  Func -> SetParameters(par);
  Func -> SetParLimits(0,1.4,1.65);
  gPlot -> Fit("Func","EM","", minx, maxx);
}
void fit_gfunctionExt(TH1F * gPlot, float minx, float maxx){
  Double_t par[5] = {1.1,0.2,0.1,1.6,0.35};
  TF1* Func = new TF1("Func", func_fitExt, 0.6,2.1,5);
  Func -> SetParameters(par);
  //Func -> SetParLimits(0,0.9,1.3);
  //Func -> SetParLimits(3,1.55,1.65);
  gPlot -> Fit("Func","EM","", minx, maxx);
}
// ====================================================================================================================================================
// polynomial background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}
// Lorentzian Peak function
Double_t gaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*(pow((x[0]-par[1])/par[2],2.0)));  
}
// Sum of background and peak function
Double_t fitSB(Double_t *x, Double_t *par) {
  return background(x,par) + gaussianPeak(x,&par[4]);
}

