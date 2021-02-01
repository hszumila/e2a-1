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

const double Me      = .000511;         // Electron     mass in GeV
const double Mp      = 0.93827;         // Proton       mass in GeV
const double Mn      = 0.93957;         // Neutron      mass in GeV
const double Mpip    = 0.13957;         // Charged pion mass in GeV
const double He3_be  = 0.00770;         // 3He binding energy in GeV
const double M_tar   = 2*Mp+Mn-He3_be;  // 3Helium  mass in GeV (target)

double max_outta_three ( double A , double B, double C );

TLorentzVector correctPn(TLorentzVector mom_n, double eIn, double eOut){
  double a,b,c;

  TVector3 mom = mom_n.Vect();
  double P = mom.Mag();
  
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    /*a = 0.06;
    b = -0.17;
    c = 0.08;*/
    a = 0.122428;
    b = -0.00330188;
    c = 0.0451578;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    /*a = 0.02;
    b = -0.05;
    c = 0.05;*/
    a = 0.173537;
    b = -0.0527054;
    c = 0.0785062;
   }
  //both
  else {
    a = 0.09;
    b = -0.25;
    c = 0.08;
    
  }

  double newMag = P - P*(a+b*P+c*P*P);
  TVector3 unit = mom.Unit();
  unit.SetMag(newMag);

  mom_n.SetVect(unit);
  
  return mom_n;
}

TLorentzVector calculatePn(double beta, TLorentzVector mom){

  TVector3 pn = mom.Vect();
  TVector3 unit = pn.Unit();
  double corrP = Mn*beta/sqrt(1.0-beta*beta);
  unit.SetMag(corrP);
  mom.SetVect(unit);
  return mom;

}



int main(){

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

   bool eppn, eppnpipi;
   const double c_cm_ns = 29.9792458;
  
  //Define Branch variables and other variables needed to determine parameters
  int nEntries;
  
  const double Mn_up   = 1.0000;		// Delta neutron mass in GeV for mass distribution cut
  const double Mn_do   = 0.8000;		// Delta neutron mass in GeV for mass distribution cut
  
  double Ebeam, M_miss_eppn, M_miss_eppnpipi, nBeta, Rcalc, cos_th_miss_ec_eppn, cos_th_miss_ec_eppnpipi;
  TVector3 EC_hit_pos, P_miss_eppn, P_miss_eppnpipi, u1;
  TLorentzVector mom_e,mom_p1,mom_p2, mom_n, E0, V4_tar, mom_pip, mom_pim;
  double sig_e_p1, sig_e_p2, sig_p1p2, sig_e_pip, sig_e_pim, sig_pip_pim, temp1, temp2,
    Z_e_p1_abs, Z_e_p2_abs, Z_p1p2_abs, Z_e_pip_abs, Z_e_pim_abs, Z_pip_pim_abs, vz_o_sig;

  //make histograms
  TH1F *h_tres = new TH1F("h_tres","EC time resolution;t_{EC}-t_{TOF} [ns]",100,-2,2);
  TH1F *h_pdiff_tot = new TH1F("h_pdiff_tot",";p_{meas}-p_{miss} [GeV/c]",50,-1,1);
  TH1F *h_pdiff_a = new TH1F("h_pdiff_a","0.5<p_{meas}<1 GeV/c;p_{meas}-p_{miss} [GeV/c]",20,-1,1);
  TH1F *h_pdiff_b = new TH1F("h_pdiff_b","1<p_{meas}<1.4 GeV/c;p_{meas}-p_{miss} [GeV/c]",20,-1,1);
  TH1F *h_pdiff_c = new TH1F("h_pdiff_c","1.4<p_{meas}<2. GeV/c;p_{meas}-p_{miss} [GeV/c]",20,-1,1);
  TH1F *nmom_a = new TH1F("nmom_a","0.5<p_{meas}<1 GeV/c",20,0,2.5);
  TH1F *nmom_b = new TH1F("nmom_b","1<p_{meas}<1.4 GeV/c",20,0,2.5);
  TH1F *nmom_c = new TH1F("nmom_c","1.4<p_{meas}<2. GeV/c",20,0,2.5); 
  TH2F *h_pdiffVp = new TH2F("h_pdiffVp",";p_{meas}-p_{miss} [GeV/c];p_{meas} [GeV/c]",50,-3,3,50,0,2.5);
  TH2F *h_pmVpm = new TH2F("h_pmVpm",";p_{meas} [GeV/c];p_{miss} [GeV/c]",50,0,3,50,0,3);
  TH1F *h_rmDot = new TH1F("h_rmDot",";cos(#theta)",100,0,180);
  TH1F *h_nmom = new TH1F("h_nmom",";p_{n,meas}",100,-5,5);
  TH1F *h_pmom = new TH1F("h_pmom",";p_{p,meas}",100,-5,5);
  TH1F *h_miss = new TH1F("h_miss",";p_{miss}",100,-5,5);
  TH1F *h_mmiss = new TH1F("h_mmiss",";m_{miss}",100,0.,1.5);

  //read in the data and make output file  
  TFile *f = new TFile("../../../data/n_eff/skim_He3_pp_2261.root");
  TFile *fout = new TFile("output/neff_He3_2261.root","RECREATE");
  /*
  TCanvas *canvas = new TCanvas("canvas","output plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/neff_He3.pdf";
  */
  TTree * t = (TTree*)f->Get("T");
  TString outputpdf = Form("out_neff/output_nMomentumRes_2261_plc");

 // --------------------------------------------------------------------------
  int nParticles, nProtons, nNeutrons, nPiplus, nPiminus;
  int Part_type [20], stat_dc  [20], stat_sc[20], stat_ec[20];
  double t0, QSq, xB, Nu, start_time;
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
  t->SetBranchAddress("start_time" , &start_time );

  Fiducial fid_params(2261,2250,5996,"3He",true);

  //-----------------------------
  //make histograms
  // Vertex difference cut
  TH1F *h1_delta_vz_e_p1    = new TH1F("h1_delta_vz_e_p1"   ,"#Delta v_{Z};e - p1;Counts"     , 50,  -5., 5. );
  TH1F *h1_delta_vz_e_p2    = new TH1F("h1_delta_vz_e_p2"   ,"#Delta v_{Z};e - p2;Counts"     , 50,  -5., 5. );
  TH1F *h1_delta_vz_p1_p2   = new TH1F("h1_delta_vz_p1_p2"  ,"#Delta v_{Z};p1 - p2;Counts"    , 50,  -5., 5. );
  TH1F *h1_delta_vz_e_pip   = new TH1F("h1_delta_vz_e_pip"  ,"#Delta v_{Z};e - pi+;Counts"    , 50,  -5., 5. );
  TH1F *h1_delta_vz_e_pim   = new TH1F("h1_delta_vz_e_pim"  ,"#Delta v_{Z};e - pi-;Counts"    , 50,  -5., 5. );
  TH1F *h1_delta_vz_pip_pim = new TH1F("h1_delta_vz_pip_pim","#Delta v_{Z};pi+ - pi-;Counts"  , 50,  -5., 5. );
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
  
  
  //pmiss vs pn both (uncorrected)
  TH2D* h_pmVpn_both = new TH2D("h_pmVpn_both","all EC;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_both = new TH2D("h_pdiffVpn_both","both EC;p_{n} [GeV/c];(p_{n}-p_{miss})/p_{n} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_both1 = new TH1F("h1_pmdiff_both1","0.5<p_{miss}<1.2 GeV (EC both, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_both2 = new TH1F("h1_pmdiff_both2","1.2<p_{miss}<2.2 GeV (EC both, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  
  //pmiss vs pn inner EC (uncorrected)
  TH2D* h_pmVpn_inner = new TH2D("h_pmVpn_inner","inner EC;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_inner = new TH2D("h_pdiffVpn_inner","inner EC;p_{n} [GeV/c];(p_{n}-p_{miss})/p_{n} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_inner1 = new TH1F("h1_pmdiff_inner1","0.5<p_{miss}<1.2 GeV (EC in, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_inner2 = new TH1F("h1_pmdiff_inner2","1.2<p_{miss}<2.2 GeV (EC in, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  //TH1F *h1_pmdiff_inner3 = new TH1F("h1_pmdiff_inner3","1.5<p_{miss}<2.1 GeV (EC in, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);

  //pmiss vs pn outer EC (uncorrected)
  TH2D* h_pmVpn_outer = new TH2D("h_pmVpn_outer","outer EC;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_outer = new TH2D("h_pdiffVpn_outer","outer EC;p_{n} [GeV/c];(p_{n}-p_{miss})/p_{n} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_outer1 = new TH1F("h1_pmdiff_outer1","0.5<p_{miss}<1.2 GeV (EC out, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_outer2 = new TH1F("h1_pmdiff_outer2","1.2<p_{miss}<2.2 GeV (EC out, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  //TH1F *h1_pmdiff_outer3 = new TH1F("h1_pmdiff_outer3","1.5<p_{miss}<2.1 GeV (EC out, uncorr);p_{miss}-p_{n} [GeV/c]",30,-1,1);

  //pmiss vs pn both (corrected)
  TH2D* h_pmVpn_bothC = new TH2D("h_pmVpn_bothC","all EC, corrected;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_bothC = new TH2D("h_pdiffVpn_bothC","both EC;p_{n} [GeV/c];(p_{n}-p_{miss})/p_{n} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_both1C = new TH1F("h1_pmdiff_both1C","0.5<p_{miss}<1.2 GeV (EC both, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_both2C = new TH1F("h1_pmdiff_both2C","1.2<p_{miss}<2.2 GeV (EC both, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
 
  //pmiss vs pn inner EC (corrected)
  TH2D* h_pmVpn_innerC = new TH2D("h_pmVpn_innerC","inner EC, corrected;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_innerC = new TH2D("h_pdiffVpn_innerC","inner EC, corr;p_{n} [GeV/c];(p_{n,corr}-p_{miss})/p_{n,corr} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_inner1C = new TH1F("h1_pmdiff_inner1C","0.5<p_{miss}<1.2 GeV (EC in, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_inner2C = new TH1F("h1_pmdiff_inner2C","1.2<p_{miss}<2.2 GeV (EC in, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  //TH1F *h1_pmdiff_inner3C = new TH1F("h1_pmdiff_inner3C","1.5<p_{miss}<2.1 GeV (EC in, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  
  //pmiss vs pn outer EC (corrected)
  TH2D* h_pmVpn_outerC = new TH2D("h_pmVpn_outerC","outer EC, corrected;p_{n} [GeV/c];p_{miss} [GeV/c]",100,0,2.3,100,0,2.3);
  TH2D* h_pdiffVpn_outerC = new TH2D("h_pdiffVpn_outerC","outer EC, corr;p_{n} [GeV/c];(p_{n,corr}-p_{miss})/p_{n,corr} [GeV/c]",100,0,2.3,100,-0.8,0.8);
  TH1F *h1_pmdiff_outer1C = new TH1F("h1_pmdiff_outer1C","0.5<p_{miss}<1.2 GeV (EC out, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  TH1F *h1_pmdiff_outer2C = new TH1F("h1_pmdiff_outer2C","1.2<p_{miss}<2.2 GeV (EC out, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);
  //TH1F *h1_pmdiff_outer3C = new TH1F("h1_pmdiff_outer3C","1.5<p_{miss}<2.1 GeV (EC out, corr);p_{miss}-p_{n} [GeV/c]",30,-1,1);

  //path length studies:
  TH1F *h1_Rec_measured = new TH1F("h_Rec_measured","R_{EC} using EC;R_{EC}",100,400,600);
  TH1F *h1_Rec_skim = new TH1F("h_Rec_skim","R_{EC,skim plc} from skim #beta;R_{EC,skim plc}",100,400,600);
  TH1F *h1_Rec_calc = new TH1F("h_Rec_calc","R_{EC,calc} using P_{miss};R_{EC,calc}",100,400,600);
  TH2D *h_RmVRc = new TH2D("h_RmVRc",";R_{EC,calc} (using pmiss) [cm];R_{EC} (measured) [cm]",100,400,600,50,500,600);
  TH2D *h_RmVRs = new TH2D("h_RmVRs",";R_{EC,skim} (using skim beta) [cm];R_{EC} (measured) [cm]",100,400,600,50,500,600);
  TH1F *h_Rdiff_inner = new TH1F("h_Rdiff_inner","inner EC;R_{calculated}-R_{measured} [cm]",100,-100,50);
  TH1F *h_Rdiff_outer = new TH1F("h_Rdiff_outer","outer EC;R_{calculated}-R_{measured} [cm]",100,-100,50);
  TH1F *h_Rdiff_both = new TH1F("h_Rdiff_both","both EC layers;R_{calculated}-R_{measured} [cm]",100,-100,50);
  TH1F *h_Rdiff_compared = new TH1F("h_Rdiff_compared","all EC, calc using skim #beta;R_{calculated}-R_{measured} [cm]",100,-100,50);
  TH1F *h_RcorrPmissDiff = new TH1F("h_RcorrPmissDiff","p_{miss}-p_{n,Rcorr}",100,-1,1);
  TH1F *h_noRcorrPmissDiff = new TH1F("h_noRcorrPmissDiff","p_{miss}-p_{n,noRcorr}",100,-1,1);
  TH1F *h_checkRcorrPmissDiff = new TH1F("h_checkRcorrPmissDiff","p_{miss}-p_{n,Rcorr}",100,-5,5);


  
  //loop events a first time to get the vertex cuts
  Ebeam = 2.2;
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
      else if(Part_type[i] == 2112){n1 = i;}

      //else if ((charge[i] == 0)&&(n1==0)){
      //	n1 = i; 
      //}
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

  //loop exclusive events ppn
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
      else if(Part_type[i] == 2112){n1 = i;}
	//((charge[i] == 0)&&(n1==0))	n1 = i;
    }
    
    if (p1==0||p2==0){cout << "No two protons in this event" << endl;	continue;}
    
    // --------------------------------------------------------------------------
    mom_e .SetXYZM(mom_x[ 0],mom_y[ 0],mom_z[ 0],Me);               // Electron
    mom_p1.SetXYZM(mom_x[p1],mom_y[p1],mom_z[p1],Mp);               // Proton 1
    mom_p2.SetXYZM(mom_x[p2],mom_y[p2],mom_z[p2],Mp);               // Proton 2

    // Calculate missing mass and momentum	
    M_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2).M   ();
    P_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2).Vect();
  
    Z_e_p1_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p1])/sig_e_p1;
    Z_e_p2_abs = TMath::Abs(vtx_z_cor[ 0]-vtx_z_cor[p2])/sig_e_p2;
    Z_p1p2_abs = TMath::Abs(vtx_z_cor[p1]-vtx_z_cor[p2])/sig_p1p2;
    vz_o_sig = max_outta_three(Z_e_p1_abs,Z_e_p2_abs,Z_p1p2_abs);

    if((M_miss_eppn < Mn_up)&&(M_miss_eppn > Mn_do)){
      if(vz_o_sig<3.){
	if(P_miss_eppn.Mag() > 0.5) {
	  if (n1!=0){ //No neutral hit in the EC
	    if( (stat_ec[n1]>0) && (stat_dc[n1]<=0) && (stat_sc[n1]<=0) ){	    
	      if(beta[n1] < 0.95){
		EC_hit_pos.SetXYZ(ec_x[n1],ec_y[n1],ec_z[n1]);
		cos_th_miss_ec_eppn = P_miss_eppn*EC_hit_pos/P_miss_eppn.Mag()/EC_hit_pos.Mag();	
		if(cos_th_miss_ec_eppn > 0.995) {
		  if((ec_tot[n1]>0)&&(ec_time[n1]>0)){

		    mom_n.SetXYZM(mom_x[n1],mom_y[n1],mom_z[n1],Mn); 
		    h_pmVpn_both->Fill(mom_n.P(),P_miss_eppn.Mag());

		    //cout<<"original mom: "<<mom_n.Vect().Mag()<<" ";
		    //mom_n = calculatePn(beta[n1],mom_n);
		    //cout<<"with plc mom: "<<mom_n.Vect().Mag()<<endl;

		    TLorentzVector mom_nCorr = correctPn(mom_n,ec_in[n1],ec_out[n1]);
		    h_pmVpn_bothC->Fill(mom_nCorr.P(),P_miss_eppn.Mag());

		    h1_Rec_measured->Fill(sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		    h1_Rec_skim->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns));
		    h1_Rec_calc->Fill(P_miss_eppn.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppn.Mag(),2.0)));		    
		    h_RmVRc->Fill(P_miss_eppn.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppn.Mag(),2.0)),sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		    h_RmVRs->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns),sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		    h_Rdiff_compared->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));


		    double dR = 0;
		    
		    
		    //inner EC
		    if (ec_in[n1]>0.01 && ec_out[n1]<0.01){
		      h_pmVpn_inner->Fill(mom_n.P(),P_miss_eppn.Mag());
		      h_pmVpn_innerC->Fill(mom_nCorr.P(),P_miss_eppn.Mag());
		      h_pdiffVpn_inner->Fill(mom_n.P(),(mom_n.P()-P_miss_eppn.Mag())/mom_n.P());
		      h_pdiffVpn_innerC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppn.Mag())/mom_nCorr.P());
		      h_Rdiff_inner->Fill(P_miss_eppn.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppn.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -44.3;
		      if (P_miss_eppn.Mag()>0.5 && P_miss_eppn.Mag()<1.2){
			h1_pmdiff_inner1->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_inner1C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppn.Mag()>=1.2 && P_miss_eppn.Mag()<2.2){
			h1_pmdiff_inner2->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_inner2C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }
		      /*else if (P_miss_eppn.Mag()>=1.5 && P_miss_eppn.Mag()<2.1){
			h1_pmdiff_inner3->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_inner3C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
			}*/
		    }
		    //outer EC
		    else if (ec_in[n1]<0.01 && ec_out[n1]>0.01){
		      h_pmVpn_outer->Fill(mom_n.P(),P_miss_eppn.Mag());
		      h_pmVpn_outerC->Fill(mom_nCorr.P(),P_miss_eppn.Mag());
		      h_pdiffVpn_outer->Fill(mom_n.P(),(mom_n.P()-P_miss_eppn.Mag())/mom_n.P());
		      h_pdiffVpn_outerC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppn.Mag())/mom_nCorr.P());
		      h_Rdiff_outer->Fill(P_miss_eppn.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppn.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -59.5;
		      if (P_miss_eppn.Mag()>0.5 && P_miss_eppn.Mag()<1.2){
			h1_pmdiff_outer1->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_outer1C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppn.Mag()>=1.2 && P_miss_eppn.Mag()<2.2){
			h1_pmdiff_outer2->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_outer2C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }
		      /* else if (P_miss_eppn.Mag()>=1.5 && P_miss_eppn.Mag()<2.1){
			h1_pmdiff_outer3->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_outer3C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
			}*/
		    }

		    //both EC
		    else if(ec_in[n1]>0.01 && ec_out[n1]>0.01){
		      h_pmVpn_both->Fill(mom_n.P(),P_miss_eppn.Mag());
		      h_pmVpn_bothC->Fill(mom_nCorr.P(),P_miss_eppn.Mag());
		      h_pdiffVpn_both->Fill(mom_n.P(),(mom_n.P()-P_miss_eppn.Mag())/mom_n.P());
		      h_pdiffVpn_bothC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppn.Mag())/mom_nCorr.P());
		      h_Rdiff_both->Fill(P_miss_eppn.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppn.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -47.2;
		      if (P_miss_eppn.Mag()>0.5 && P_miss_eppn.Mag()<1.2){
			h1_pmdiff_both1->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_both1C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppn.Mag()>=1.2 && P_miss_eppn.Mag()<2.2){
			h1_pmdiff_both2->Fill(P_miss_eppn.Mag()-mom_n.P());
			h1_pmdiff_both2C->Fill(P_miss_eppn.Mag()-mom_nCorr.P());
		      }

		    }

		    double betaNew = (sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)) + dR)/((ec_time[n1]-start_time)*c_cm_ns);
		    double pnMag = Mn*betaNew/sqrt(1.-betaNew*betaNew);
		    h_RcorrPmissDiff->Fill(P_miss_eppn.Mag()-pnMag);
		    h_checkRcorrPmissDiff->Fill(P_miss_eppn.Mag()-pnMag);
		    h_noRcorrPmissDiff->Fill(P_miss_eppn.Mag()-mom_n.P());
		  }}}}}}}}}
  
		    
  //loop inclusive event pppi+pi-n
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
      //      else if ((charge[i] == 0)&&(n1==0))	n1 = i;
      else if(Part_type[i] == 2112){n1 = i;}

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
    
    Z_e_p1_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p1])/sig_e_p1   ;
    Z_e_p2_abs    = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[ p2])/sig_e_p2   ;
    Z_p1p2_abs    = TMath::Abs(vtx_z_cor[ p1]-vtx_z_cor[ p2])/sig_p1p2   ;
    Z_e_pip_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pip])/sig_e_pip  ;
    Z_e_pim_abs   = TMath::Abs(vtx_z_cor[  0]-vtx_z_cor[pim])/sig_e_pim  ;
    Z_pip_pim_abs = TMath::Abs(vtx_z_cor[pip]-vtx_z_cor[pim])/sig_pip_pim;
    temp1    = max_outta_three(Z_e_p1_abs ,Z_e_p2_abs ,Z_p1p2_abs   );
    temp2    = max_outta_three(Z_e_pip_abs,Z_e_pim_abs,Z_pip_pim_abs);
    vz_o_sig = TMath::Max(temp1,temp2);

    if((M_miss_eppnpipi < Mn_up)&&(M_miss_eppnpipi > Mn_do)){
      if(vz_o_sig<3.){
	if(P_miss_eppnpipi.Mag() > 0.5){
	  if (n1!=0){ //No neutral hit in the EC
	    
	    if( (stat_ec[n1]>0) && (stat_dc[n1]<=0) && (stat_sc[n1]<=0) ){
	      if(beta[n1] < 0.95){
		EC_hit_pos.SetXYZ(ec_x[n1],ec_y[n1],ec_z[n1]);
		cos_th_miss_ec_eppnpipi = P_miss_eppnpipi*EC_hit_pos/P_miss_eppnpipi.Mag()/EC_hit_pos.Mag();
		if(cos_th_miss_ec_eppnpipi > 0.995) {
		  if((ec_tot[n1]>0)&&(ec_time[n1]>0)){

		    mom_n.SetXYZM(mom_x[n1],mom_y[n1],mom_z[n1],Mn); 
		    h_pmVpn_both->Fill(mom_n.P(),P_miss_eppnpipi.Mag());

		    //mom_n = calculatePn(beta[n1],mom_n);
		    
		    TLorentzVector mom_nCorr = correctPn(mom_n,ec_in[n1],ec_out[n1]);
		    h_pmVpn_bothC->Fill(mom_nCorr.P(),P_miss_eppnpipi.Mag());

		    h1_Rec_measured->Fill(sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		    h1_Rec_skim->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns));
		    h1_Rec_calc->Fill(P_miss_eppnpipi.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppnpipi.Mag(),2.0)));
		    
		    h_RmVRc->Fill(P_miss_eppnpipi.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppnpipi.Mag(),2.0)),sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		    h_RmVRs->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns),sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));		   
		   
		    h_Rdiff_compared->Fill(beta[n1]*((ec_time[n1]-start_time)*c_cm_ns) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));

		    double dR = 0;
		    
		    //inner EC
		    if (ec_in[n1]>0.01 && ec_out[n1]<0.01){
		      h_pmVpn_inner->Fill(mom_n.P(),P_miss_eppnpipi.Mag());
		      h_pmVpn_innerC->Fill(mom_nCorr.P(),P_miss_eppnpipi.Mag());
		      h_pdiffVpn_inner->Fill(mom_n.P(),(mom_n.P()-P_miss_eppnpipi.Mag())/mom_n.P());
		      h_pdiffVpn_innerC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppnpipi.Mag())/mom_nCorr.P());
		      h_Rdiff_inner->Fill(P_miss_eppnpipi.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppnpipi.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -44.3;
		      if (P_miss_eppnpipi.Mag()>0.5 && P_miss_eppnpipi.Mag()<1.2){
			h1_pmdiff_inner1->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_inner1C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppnpipi.Mag()>=1.2 && P_miss_eppnpipi.Mag()<2.2){
			h1_pmdiff_inner2->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_inner2C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }
		      /*else if (P_miss_eppnpipi.Mag()>=1.5 && P_miss_eppnpipi.Mag()<2.1){
			h1_pmdiff_inner3->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_inner3C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
			}*/
		    }
		    //outer EC
		    else if (ec_in[n1]<0.01 && ec_out[n1]>0.01){
		      h_pmVpn_outer->Fill(mom_n.P(),P_miss_eppnpipi.Mag());
		      h_pmVpn_outerC->Fill(mom_nCorr.P(),P_miss_eppnpipi.Mag());
		      h_pdiffVpn_outer->Fill(mom_n.P(),(mom_n.P()-P_miss_eppnpipi.Mag())/mom_n.P());
		      h_pdiffVpn_outerC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppnpipi.Mag())/mom_nCorr.P());
		      h_Rdiff_outer->Fill(P_miss_eppnpipi.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppnpipi.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -59.5;

		      if (P_miss_eppnpipi.Mag()>0.5 && P_miss_eppnpipi.Mag()<1.2){
			h1_pmdiff_outer1->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_outer1C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppnpipi.Mag()>=1.2 && P_miss_eppnpipi.Mag()<2.2){
			h1_pmdiff_outer2->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_outer2C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }
		      /*else if (P_miss_eppnpipi.Mag()>=1.5 && P_miss_eppnpipi.Mag()<2.1){
			h1_pmdiff_outer3->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_outer3C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
			}*/
		    }

		    //both EC
		    else if(ec_in[n1]>0.01 && ec_out[n1]>0.01){
		      h_pmVpn_both->Fill(mom_n.P(),P_miss_eppnpipi.Mag());
		      h_pmVpn_bothC->Fill(mom_nCorr.P(),P_miss_eppnpipi.Mag());
		      h_pdiffVpn_both->Fill(mom_n.P(),(mom_n.P()-P_miss_eppnpipi.Mag())/mom_n.P());
		      h_pdiffVpn_bothC->Fill(mom_n.P(),(mom_nCorr.P()-P_miss_eppnpipi.Mag())/mom_nCorr.P());
		      h_Rdiff_both->Fill(P_miss_eppnpipi.Mag()*c_cm_ns*(ec_time[n1]-start_time)/sqrt(pow(Mn,2.)+pow(P_miss_eppnpipi.Mag(),2.0)) - sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)));
		      dR = -47.2;

		      if (P_miss_eppnpipi.Mag()>0.5 && P_miss_eppnpipi.Mag()<1.2){
			h1_pmdiff_both1->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_both1C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }
		      else if (P_miss_eppnpipi.Mag()>=1.2 && P_miss_eppnpipi.Mag()<2.2){
			h1_pmdiff_both2->Fill(P_miss_eppnpipi.Mag()-mom_n.P());
			h1_pmdiff_both2C->Fill(P_miss_eppnpipi.Mag()-mom_nCorr.P());
		      }

		    }
		    double betaNew = (sqrt(pow(ec_x[n1],2.)+pow(ec_y[n1],2.)+pow(ec_z[n1]-vtx_z_cor[0],2.0)) + dR)/((ec_time[n1]-start_time)*c_cm_ns);
		    double pnMag = Mn*betaNew/sqrt(1.-betaNew*betaNew);
		    h_RcorrPmissDiff->Fill(P_miss_eppnpipi.Mag()-pnMag);
		    h_checkRcorrPmissDiff->Fill(P_miss_eppnpipi.Mag()-pnMag);
		    h_noRcorrPmissDiff->Fill(P_miss_eppnpipi.Mag()-mom_n.P());

		  }}}}}}}}}


  //slice and dice
  int nbins = 8;
  TH1F *h_inslice[nbins];
  TH1F *h_outslice[nbins];
  double pnslice[nbins];
  double err_pnslice[nbins];
  double mean_inner[nbins];
  double mean_outer[nbins];
  TF1 *c_in[nbins];
  TF1 *c_out[nbins];
  double err_inner[nbins];
  double err_outer[nbins];
  for (int kk=0;kk<nbins;kk++){
    pnslice[kk] = (double) 0.6 + 0.2*kk;
    cout<<"pn slice: "<<kk<<" "<<pnslice[kk]<<endl;;
    err_pnslice[kk] = 0.1;
    h_inslice[kk] = new TH1F(Form("h_inslice_%d",kk),"",30,-1,1);
    h_outslice[kk] = new TH1F(Form("h_outslice_%d",kk),"",30,-1,1);
  }

  for (int ii=0; ii<nbins; ii++){
    int lowBin = h_pdiffVpn_inner->GetXaxis()->FindBin(pnslice[ii]-0.1);
    int highBin = h_pdiffVpn_inner->GetXaxis()->FindBin(pnslice[ii]+0.1);
    h_inslice[ii] = (TH1F*)h_pdiffVpn_inner->ProjectionY(Form("sliceIn_%d",ii),lowBin,highBin);
    h_outslice[ii] = (TH1F*)h_pdiffVpn_outer->ProjectionY(Form("sliceOut_%d",ii),lowBin,highBin);
    c_in[ii] = new TF1(Form("cin%d",ii),"gaus",-1,1);
    c_out[ii] = new TF1(Form("cout%d",ii),"gaus",-1,1);
    h_inslice[ii]->Fit(c_in[ii],"QS");
    h_outslice[ii]->Fit(c_out[ii],"QS");
    mean_inner[ii] = c_in[ii]->GetParameter(1);
    mean_outer[ii] = c_out[ii]->GetParameter(1);
    err_inner[ii] = 1./sqrt(h_inslice[ii]->GetEntries());
    err_outer[ii] = 1./sqrt(h_outslice[ii]->GetEntries());
    cout<<"means: "<<ii<<" "<<mean_inner[ii]<<" "<<mean_outer[ii]<<" "<<h_inslice[ii]->GetEntries()<<" "<<h_outslice[ii]->GetEntries()<<endl;;

  }
  
 
 
  TCanvas *c1 = new TCanvas("c1","",1200,900);
  c1->Divide(2,1);
  c1->cd(1);
  h_pmVpn_inner->Draw("colz");

  c1->cd(2);
  h_pmVpn_outer->Draw("colz");
  c1->Print(outputpdf+".pdf(");

  TCanvas *c2 = new TCanvas("c2","",1200,900);
  c2->Divide(2,1);

  c2->cd(1);
  h_pmVpn_innerC->Draw("colz");

  c2->cd(2);
  h_pmVpn_outerC->Draw("colz");
  c2->Print(outputpdf+".pdf(");

  TCanvas *c22 = new TCanvas("c22","",1200,900);
  c22->Divide(2,1);
  c22->cd(1);
  h_pmVpn_both->Draw("colz");
  c22->cd(2);
  h_pmVpn_bothC->Draw("colz");
  c22->Print(outputpdf+".pdf(");
  
  TCanvas *c3 = new TCanvas("c3","",1200,900);
  c3->Divide(2,1);
  c3->cd(1);
  h1_pmdiff_inner1->Draw("");
  c3->cd(2);
  h1_pmdiff_inner1C->Draw("");
  c3->Print(outputpdf+".pdf(");

  TCanvas *c4 = new TCanvas("c4","",1200,900);
  c4->Divide(2,1);
  c4->cd(1);
  h1_pmdiff_inner2->Draw("");
  c4->cd(2);
  h1_pmdiff_inner2C->Draw("");
  c4->Print(outputpdf+".pdf(");
  /*
  TCanvas *c5 = new TCanvas("c5","",1200,900);
  c5->Divide(2,1);
  c5->cd(1);
  h1_pmdiff_inner3->Draw("");
  c5->cd(2);
  h1_pmdiff_inner3C->Draw("");
  c5->Print(outputpdf+".pdf(");
  */
  TCanvas *c6 = new TCanvas("c6","",1200,900);
  c6->Divide(2,1);
  c6->cd(1);
  h1_pmdiff_outer1->Draw("");
  c6->cd(2);
  h1_pmdiff_outer1C->Draw("");
  c6->Print(outputpdf+".pdf(");

  TCanvas *c7 = new TCanvas("c7","",1200,900);
  c7->Divide(2,1);
  c7->cd(1);
  h1_pmdiff_outer2->Draw("");
  c7->cd(2);
  h1_pmdiff_outer2C->Draw("");
  c7->Print(outputpdf+".pdf(");

   TCanvas *c66 = new TCanvas("c66","",1200,900);
  c66->Divide(2,1);
  c66->cd(1);
  h1_pmdiff_both1->Draw("");
  c66->cd(2);
  h1_pmdiff_both1C->Draw("");
  c66->Print(outputpdf+".pdf(");

  TCanvas *c77 = new TCanvas("c77","",1200,900);
  c77->Divide(2,1);
  c77->cd(1);
  h1_pmdiff_both2->Draw("");
  c77->cd(2);
  h1_pmdiff_both2C->Draw("");
  c77->Print(outputpdf+".pdf(");
  
  /*
  TCanvas *c8 = new TCanvas("c8","",1200,900);
  c8->Divide(2,1);
  c8->cd(1);
  h1_pmdiff_outer3->Draw("");
  c8->cd(2);
  h1_pmdiff_outer3C->Draw("");
  c8->Print(outputpdf+".pdf)");
  */

  TCanvas *c9 = new TCanvas("c9","",1200,900);
  c9->Divide(2,1);
  c9->cd(1);
  h_pdiffVpn_inner->Draw("colz");
  c9->cd(2);
  h_pdiffVpn_innerC->Draw("colz");
  c9->Print(outputpdf+".pdf(");

  TCanvas *c10 = new TCanvas("c10","",1200,900);
  c10->Divide(2,1);
  c10->cd(1);
  h_pdiffVpn_both->Draw("colz");
  c10->cd(2);
  h_pdiffVpn_bothC->Draw("colz");
  c10->Print(outputpdf+".pdf(");

  TCanvas *c11 = new TCanvas("c11","",1200,900);
  c11->Divide(2,1);
  c11->cd(1);
  h_pdiffVpn_outer->Draw("colz");
  c11->cd(2);
  h_pdiffVpn_outerC->Draw("colz");
  c11->Print(outputpdf+".pdf(");

  TCanvas *cnew[nbins];
  for (int ii=0;ii<nbins;ii++){
    cnew[ii] = new TCanvas(Form("cnew_%d",ii),"",1200,900);
    cnew[ii]->Divide(2,1);
    cnew[ii]->cd(1);
    h_inslice[ii]->Draw();

    cnew[ii]->cd(2);
    h_outslice[ii]->Draw();
    cnew[ii]->Print(outputpdf+".pdf(");
  }
  //cnew[0]->Print(outputpdf+".pdf)");
  
  TGraphErrors *gin = new TGraphErrors(nbins,pnslice,mean_inner,err_pnslice,err_inner);
  TGraphErrors *gout = new TGraphErrors(nbins,pnslice,mean_outer,err_pnslice,err_outer);

  TCanvas *cc = new TCanvas("cc","",1200,900);
  cc->cd();
  gin->Draw("AP");
  cc->Update();
  TF1 *innerfit = new TF1("innerfit","[0]+[1]*x+[2]*x*x",0.65,2.1);
  gin->Fit(innerfit,"R");
  cc->SaveAs("gin.C");

  TCanvas *cc1 = new TCanvas("cc1","",1200,900);
  cc1->cd();
  gout->Draw("AP");
  cc1->Update();
  TF1 *outerfit = new TF1("outerfit","[0]+[1]*x+[2]*x*x",0.65,2.1);
  gout->Fit(outerfit,"R");
  cc1->SaveAs("gout.C");

  TCanvas *cr1 = new TCanvas("cr1","",1200,900);
  cr1->cd();
  h1_Rec_measured->Draw();
  cr1->Print(outputpdf+".pdf(");

  TCanvas *cr2 = new TCanvas("cr2","",1200,900);
  cr2->cd();
  h1_Rec_skim->Draw();
  cr2->Print(outputpdf+".pdf(");

  TCanvas *cr3 = new TCanvas("cr3","",1200,900);
  cr3->cd();
  h1_Rec_calc->Draw();
  cr3->Print(outputpdf+".pdf(");

  TCanvas *cr4 = new TCanvas("cr4","",1200,900);
  cr4->cd();
  h_RmVRc->Draw("colz");
  cr4->Print(outputpdf+".pdf(");

  TCanvas *cr5 = new TCanvas("cr5","",1200,900);
  cr5->cd();
  h_RmVRs->Draw("colz");
  cr5->Print(outputpdf+".pdf(");

  TCanvas *cr6 = new TCanvas("cr6","",1200,900);
  cr6->cd();
  TF1 *fr_inner = new  TF1("fr_inner","gaus",-80,0);
  h_Rdiff_inner->Fit(fr_inner,"R");
  h_Rdiff_inner->Draw();  
  cr6->Print(outputpdf+".pdf(");


  TCanvas *cr7 = new TCanvas("cr7","",1200,900);
  cr7->cd();
  TF1 *fr_outer = new  TF1("fr_outer","gaus",-90,-30);
  h_Rdiff_outer->Fit(fr_outer,"R");
  h_Rdiff_outer->Draw();
  cr7->Print(outputpdf+".pdf(");

  TCanvas *cr8 = new TCanvas("cr8","",1200,900);
  cr8->cd();
  TF1 *fr_both = new  TF1("fr_both","gaus",-80,0);
  h_Rdiff_both->Fit(fr_both,"R");
  h_Rdiff_both->Draw();
  cr8->Print(outputpdf+".pdf(");

  TCanvas *cr9 = new TCanvas("cr9","",1200,900);
  cr9->cd();
  h_Rdiff_compared->Draw();
  cr9->Print(outputpdf+".pdf(");


  TCanvas *cr10 = new TCanvas("cr10","",1200,900);
  cr10->cd();
  h_noRcorrPmissDiff->Draw();
  cr10->Print(outputpdf+".pdf(");

  TCanvas *cr11 = new TCanvas("cr11","",1200,900);
  cr11->cd();
  h_RcorrPmissDiff->Draw();
  cr11->Print(outputpdf+".pdf(");

  TCanvas *cr12 = new TCanvas("cr12","",1200,900);
  cr12->cd();
  h_checkRcorrPmissDiff->Draw();
  cr12->Print(outputpdf+".pdf)");

}
double max_outta_three ( double A , double B, double C ){
  
  double temp_max, max_val;
  
  temp_max = TMath::Max(A,B       );
  max_val  = TMath::Max(C,temp_max);
  
  return max_val;
}
