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

#include "Fiducial.h"


//define some constants
const double amu     = 0.931494;
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
const double cc      = 2.99792458E+8;

//SRC cuts defined:
const double vtx_min_C = 4.0;//cm
const double vtx_max_C = 7.0;//cm
const double vtx_min_He = -2.5;//cm
const double vtx_max_He = -0.5;//cm
const double xb_cut = 1.1;
const double thetapq_cut = 25.0;//deg
const double pq_min = 0.62;
const double pq_max_n = 1.1;
const double pq_max_p = 0.96;
const double mmiss_nomCut = 1.1;
const double pmiss_nomCut = 0.3;
//check if using Q2 cut!!! As is the case for 2.2 GeV data

double fn_Mmiss(double omega, double Mnuc, double Enuc, TVector3 q_vector, TVector3 p_vector){
  return sqrt(pow(omega+2*Mnuc-Enuc,2.0) - (q_vector - p_vector).Mag2());
}

double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;								// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;								// Missing energy
}

double f_delta_p(double p){
  // Returns value of delta_p as a function of p.
  // This function is used for neutron momentum resolution.
  // 0.392319c	-0.0967682	0.022168
  
  //4.4, 2.2 GeV 3He:
  double ec_extr_res = 0.098159;//0.091479;//-0.09;//0967682;
  double dt = 0.3532;//0.3292;//0.392319;//0.625                 //ec_time_res;
  const double d      = 516.9152;       //cm, distance from interaction vtx to EC
  const double Mn     = 0.939565378;    // Neutron mass in GeV
  const double c_cm_s = 2.99792458E+10; // cm/s, speed of light
  const double ns_2_s = 1E-09;
  
  double val = p*p/(Mn*Mn)*sqrt(Mn*Mn+p*p)*dt/d*c_cm_s*ns_2_s;
  return sqrt(val*val/p/p+ec_extr_res*ec_extr_res)*p;
}

double correctPn(double P, double eIn, double eOut){
  double a,b,c;
  
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    a = 0.06;
    b = -0.17;
    c = 0.08;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    a = 0.02;
    b = -0.05;
    c = 0.05;
  }
  //both
  else {
    a = 0.09;
    b = -0.25;
    c = 0.08;
  }

  P = P - P*(a+b*P+c*P*P);
  return P;
}

double f_delta_p2(double P, double eIn, double eOut){
  double a,b;
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    a = 0.03;
    b = 0.06;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    a = 0.03;
    b = 0.05;
  }
  //both
  else {
    a = 0.03;
    b = 0.04;
  }
  
  return P*(a+b*P); 
}

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

TVector3 P_Coulomb_Corrected(TVector3 P_Uncorrected, int Z, int A){
	float Correction_GeV = Coulomb_Corr_MeV( Z , A )/1000.;
	double E_uncorr = sqrt(P_Uncorrected.Mag2()+pow(Mp,2.0));
	float E_corr  = E_uncorr - Correction_GeV;
	float New_Mag = sqrt(E_corr*E_corr-Mp*Mp);

	TVector3 unit = P_Uncorrected.Unit();
	unit.SetMag(New_Mag);

	return unit;
}

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

//input is 3He, 4He, or 12C
int main(int argc, char ** argv){

  //////////////////////////////////
  // Set some target dep parameters
  /////////////////////////////////
  TString target = argv[1];
  cout<<"Target is: "<<target<<endl;

  double tgt_min, tgt_max, tgt_M;
  int tgt_Z, tgt_A;
  if (target=="12C"){tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C;}
  else if (target=="3He"){tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 3; tgt_M = M3He;}
  else if (target=="4He"){tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 4; tgt_M = M4He;}

  cout<<"target min: "<<tgt_min<<" target max: "<<tgt_max<<endl;
  //Change for fiducial cuts
  Fiducial fid_params(4461,2250,5996,"4He",true);

  //read in the data and make output file  
  //TFile *f = new TFile("../../../data/19jun20/skim_C12_p_Jun20.root");
  TFile *f = new TFile("../../../data/19jun20/skim_He4_p_Jun20.root");
  TFile *fout = new TFile("output/out_He4_p.root","RECREATE");

  std::string pdf_file_name = "output/out_He4_p_Msm.pdf";

  
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  
  //////////////////////////////////
  // Plots/////////////////////////
  /////////////////////////////////
  //plot theta vs phi before
  TH2F *hfid_pre = new TH2F("hfid_pre","no p fid cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot theta vs phi after
  TH2F *hfid_post= new TH2F("hfid_post","with p fid cut & EC 10cm cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot EC x vs y before
  TH2F *hec_pre = new TH2F("hec_pre","no fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
  //plot EC x vs y after
  TH2F *hec_post = new TH2F("hec_post","10cm fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
   //plot Q^2
  TH1F* h_q2 =new TH1F("h_q2",";Q^{2} (with lead cuts);",100,0,6);
  //plot xB
  TH1F* h_xb =new TH1F("h_xb",";x_{B} (with lead cuts);",100,0,2);
  //plot nu
  TH1F* h_nu =new TH1F("h_nu",";#nu [GeV] (with lead cuts);",100,0,3);
  //plot e- theta
  TH1F* h_etheta =new TH1F("h_etheta",";#theta_{e-} [deg] (with lead cuts);",100,0,45);
  //plot e- momentum
  TH1F* h_emom =new TH1F("h_emom",";P_{e-} [GeV/c] (with lead cuts);",100,0,6); 
  //plot pmiss and xb
  TH2D * hist_pmiss_xb = new TH2D("pmiss_xb","Events before proton cuts;pmiss;Xb;Counts",28,0.3,1.0,40,1,3);
  //plot theta_pq and p/q
  TH2D * hist_pq_thetapq = new TH2D("hist_pq_thetapq","Events before proton cuts;p/q;#theta_{pq} [deg];Counts",100,0,1.5,100,0,70);
  //plot theta_pq and p/q
  TH2D * hist_pq_thetapq_smear = new TH2D("hist_pq_thetapq_smear","Smeared events before lead cuts;p_{smeared}/q;#theta_{pq} [deg];Counts",100,0,1.5,100,0,70);
  //plot z vertex
  TH1D * h_zvtx_precut = new TH1D("h_zvtx_precut","Corrected z vertex with cut;z vertex [cm]",100,-10,10);
  //plot the missing mass
  TH1D * h_mmiss = new TH1D("h_mmiss","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  TH1D * h_mmissSmear = new TH1D("h_mmissSmear","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  TH1D * h_mmissCut = new TH1D("h_mmissCut","Missing mass with lead cuts;missing mass [GeV/c^{2}]",100,0,3);
  //plot the missing energy
  TH1D * h_emiss = new TH1D("h_emiss","Missing energy;missing energy",100,-1,2);
  TH1D * h_emissSmear = new TH1D("h_emissSmear","Missing energy;missing energy",100,-0.3,2);
  TH1D * h_emissCut = new TH1D("h_emissCut","Missing energy with lead cuts;missing energy",100,-1,2);
  //plot pmiss
  TH1F* h_pmiss =new TH1F("h_pmiss",";P_{miss} [GeV/c];",100,0,1.5);
  TH1F* h_pmissSmear =new TH1F("h_pmissSmear",";P_{miss} [GeV/c];",100,0,1.5);
  TH1F* h_pmissCut =new TH1F("h_pmissCut","Pmiss with normal SRC cuts;P_{miss} [GeV/c] (eg2 in pink);",100,0,1.5);
  //plot nuclean momentum
  TH1F* h_nmom =new TH1F("h_nmom",";P_{lead N} [GeV/c];",100,0,3);
  TH1F* h_nmomSmear =new TH1F("h_nmomSmear",";P_{lead N} [GeV/c];",100,0,3);
  TH1F* h_nmomCut =new TH1F("h_nmomCut","P_{lead N} [GeV/c] with lead cuts;P_{lead N} [GeV/c];",100,0,3);
  //plot theta_n
  TH1F* h_ntheta =new TH1F("h_ntheta",";#theta_{lead N} [deg];",100,0,90);
  TH1F* h_nthetaCut =new TH1F("h_nthetaCut","#theta_{lead N} [deg] with lead cuts;#theta_{lead N} [deg];",100,0,90);
  //plot theta_pq
  TH1F* h_thetapq =new TH1F("h_thetapq",";#theta_{lead N/q} [deg];",100,0,25);
  TH1F* h_thetapqCut =new TH1F("h_thetapqCut","#theta_{lead N/q} [deg] with lead cuts;#theta_{lead N/q} [deg];",100,0,25);
  //plot proton smearing
  TH1F* h_smear =new TH1F("h_smear",";p_{p} - p_{p,smeared} [GeV/c];",100,-1,1);
  TH1F* h_smearCut =new TH1F("h_smearCut","with lead cuts;p_{p} - p_{p,smeared} [GeV/c];",100,-1,1);
  //plot nominal quantities
  TH2F* h_emVpm =new TH2F("h_emVpm",";pmiss;emiss",100,0,0.8,100,-0.5,0.5);
  //plot smeared quantities
  TH2F* h_emVpmSm =new TH2F("h_emVpmSm","smeared;pmiss;emiss",100,0,0.8,100,-0.5,0.5);
  //plot EC time resolution
  TH1F* h_timeres = new TH1F("h_timeres","EC time resolution;t_{EC}-t_{TOF} [ns]",100,-5,5);
  //plot pmeasured - pmiss
  TH1F* h_pmMpm = new TH1F("h_pmMpm",";p_{measured}-p_{miss} [GeV]",100,-5,5);
  //plot p vs psmear
  TH2F *pvps = new TH2F("h_pvps",";p_{miss} [GeV];p_{miss, from smeared}[GeV]",100,0,1,100,0,1);

  //plot the mmiss vs pmiss:
  TH2F *h_mmvpm = new TH2F("h_mmvpm",";p_{miss} [GeV];m_{miss}[GeV]",100,0,1,100,0,1.6);
  TH2F *h_mmvpmS = new TH2F("h_mmvpmS","smeared;p_{miss} [GeV];m_{miss}[GeV]",100,0,1,100,0,1.6);

  //read in the Pmiss plots from eg2
  TFile *f2 = new TFile("../../../pmiss_eg2.root");
  TH1F *eg2_pmiss;
  TCanvas *c1 = (TCanvas*)f2->Get("Canvas_1");
  eg2_pmiss = (TH1F*)c1->GetPrimitive("h1");

  TCanvas *canvas = new TCanvas("canvas","output p plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();

  TTree * intree = (TTree*)f->Get("T");
  const int maxParticles=50;
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles],ec_in[nParticles], ec_out[nParticles];
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

  const int counts = 100000;
  double pmissO[counts];
  double pmissSm[counts];
  double mmissO[counts];
  double mmissSm[counts];
  int counter = 0;
  gRandom = new TRandom3();
  double p_over_q_O[counts];
  double p_over_q_sm[counts];
  double thetaPQ[counts];
  double thetaPQsm[counts];
  double p_over_q_prev = 0;
  double theta_prev = 0;

  for (int event =0 ; event < intree->GetEntries() ; event++)
  //for (int event =0 ; event < 100000 ; event++)
    {
      if (event%10000==0){cout << "Working on event " << event <<endl;}

      intree->GetEvent(event);

       // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //cut on xB
      if (xB<xb_cut){continue;}

      //for 2.2 gev data, cut on Q2:
      //if (QSq<1.6){continue;}

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at protons
	  if (type[part] != 2212)
	    continue;

	  //vertex cut-target dep.
	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Ratio of p/q
	  double p_over_q = this_mom.Mag() / q.Mag();
	  //if ((p_over_q < pq_min)||(p_over_q > pq_max_p)) {continue;}

	  // Angle between p and q
	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
	  if (theta_pq> thetapq_cut){continue;}

	  //if (leadingID != -1){cout<<"already lead particle!"<<endl;}
	  
	  //first proton of interest in event
	  if (leadingID==-1){p_over_q_prev=p_over_q; theta_prev=theta_pq; leadingID=part;}

	  /////////////////////////////////////////////////////////////////////
	  //Here I selected the proton most likely to pass lead if more than 1:
	  //current p/q is larger and angle is smaller than previous
	  if ((p_over_q>p_over_q_prev) && (theta_pq<theta_prev)){leadingID=part;}

	  //one passes p/q cut
	  else if ((p_over_q<pq_min || p_over_q_prev<pq_min) && (p_over_q>p_over_q_prev)){leadingID=part;}

	  //both pass p/q cut or neither pass p/q cut
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapq_cut)/thetapq_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapq_cut)/thetapq_cut)
	    {leadingID=part;}
	  	  
	  /////////////////////////////////////////////////////////////////////	  
	}

      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  	  
	  //smear the proton momentum resolution
	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p2(mom_lead.Mag(),ec_in[leadingID],ec_out[leadingID]));
	  TVector3 u1 = mom_lead.Unit();
	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

	  //recalculate the smeared factors
	  TVector3 mom_missSmear = mom_smear - q;
	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 

	  //cut at Pmiss>1:
	  //if (mom_miss.Mag()>1.0 || mom_missSmear.Mag()>1.0){continue;}
	  
	  //apply proton fiducial cuts:
	  hfid_pre->Fill(phi, theta);
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}
	  hec_pre->Fill(ec_x[leadingID],ec_y[leadingID]);

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  if(mom_smear.Mag()>2.34) continue;

	  hec_post->Fill(ec_x[leadingID],ec_y[leadingID]);	  
	  hfid_post->Fill(phi, theta);
	  h_zvtx_precut->Fill(vtxZ[leadingID]);

	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();

	  //plot the (e,e'p) selected events and count them
	  if (p_over_q>pq_min && p_over_q<pq_max_p && mmiss<mmiss_nomCut && xB>1.2 && mom_miss.Mag()<1){

	    h_pmissCut->Fill(mom_miss.Mag());

	    if (mom_miss.Mag()>pmiss_nomCut){
	      h_xb->Fill(xB);
	      h_nu->Fill(Nu);
	      h_q2->Fill(QSq);
	      h_etheta->Fill(mom_e.Theta()*180./M_PI);
	      h_emom->Fill(mom_e.Mag());
	      h_mmissCut->Fill(mmiss);
	      h_nmomCut->Fill(mom_lead.Mag());
	      h_nthetaCut->Fill(theta);
	      h_thetapqCut->Fill(mom_lead.Angle(q)*180./M_PI);
	      h_emissCut->Fill(emiss);
	      h_smearCut->Fill(mom_lead.Mag()-mom_smear.Mag());

	    }
	  }
	  //plot the looser (e,e'p) events that we use to smear
	  h_pmMpm->Fill(mom_lead.Mag()-mom_miss.Mag());
	  hist_pmiss_xb->Fill(mom_miss.Mag(),xB);
	  h_mmiss->Fill(mmiss);
	  h_pmiss->Fill(mom_miss.Mag());
	  h_nmom->Fill(mom_lead.Mag());
	  h_ntheta->Fill(theta);
	  h_thetapq->Fill(mom_lead.Angle(q)*180./M_PI);
	  h_emiss->Fill(emiss);
	  hist_pq_thetapq->Fill(mom_lead.Mag()/q.Mag(),mom_lead.Angle(q)*180./M_PI);
	  h_emVpm->Fill(mom_miss.Mag(),emiss);
	  h_timeres->Fill(sc_time[leadingID]-ec_time[leadingID]-(sc_path[leadingID]-ec_path[leadingID])/cc);

	  //plot the smeared (e,e'p) events that look like neutrons
	  pmissO[counter] = mom_miss.Mag();
	  pmissSm[counter] = mom_missSmear.Mag();
	  mmissO[counter] = mmiss;
	  mmissSm[counter] = mmissSmear;
	  p_over_q_O[counter] = mom_lead.Mag() / q.Mag();
	  p_over_q_sm[counter] = mom_smear.Mag() / q.Mag();
	  thetaPQsm[counter] = mom_smear.Angle(q);
	  thetaPQ[counter] = mom_lead.Angle(q);
	  
	  hist_pq_thetapq_smear->Fill(mom_smear.Mag()/q.Mag(),mom_smear.Angle(q)*180./M_PI);
	  h_smear->Fill(mom_lead.Mag()-mom_smear.Mag());
	  h_pmissSmear->Fill(mom_missSmear.Mag());
	  h_emissSmear->Fill(emissSmear);
	  h_mmissSmear->Fill(mmissSmear);
	  h_nmomSmear->Fill(mom_smear.Mag());
	  h_emVpmSm->Fill(mom_missSmear.Mag(),emissSmear);
	  pvps->Fill(mom_miss.Mag(),mom_missSmear.Mag());
	  if (mom_miss.Mag()>0.3&&mmiss<1.1){
	    h_mmvpmS->Fill(mom_missSmear.Mag(),mmissSmear);
	    h_mmvpm->Fill(mom_miss.Mag(),mmiss);
	  }
	  counter++;	    
	  
	}//end leading
    }//end event

  //Calculate the false positives and false negatives.
  const int n_pmiss_cut = 30;
  const int n_mmiss_cut = 7;
  //make the Pmiss cuts:
  double pmiss_cut[n_pmiss_cut];
  for (int ii=0; ii<n_pmiss_cut; ii++){
    pmiss_cut[ii] = (double)ii*2.0/100.0;
  }

  //make Mmiss cuts:
  double mmiss_val[n_mmiss_cut] = {1.1,1.15,1.175,1.2,1.225,1.25,1.3};
  double falsePos[n_mmiss_cut][n_pmiss_cut];
  double falseNeg[n_mmiss_cut][n_pmiss_cut];
  double minimum = 9999;
  int mmissCut = 0;
  int pmissCut = 0;
  double falsePosO[n_mmiss_cut][n_pmiss_cut];
  double falseNegO[n_mmiss_cut][n_pmiss_cut];
  double falsePosErr[n_mmiss_cut][n_pmiss_cut];
  double falseNegErr[n_mmiss_cut][n_pmiss_cut];

  cout<<"pmiss: "<<" mmiss: "<<" false neg: "<<" false pos: "<<" false neg error: "<<" false pos error: "<<" fneg-fpos: "<<" err(fneg-fpos): "<<endl;
  
  //loop missing mass
  for (int ii=0;ii<n_mmiss_cut; ii++){
    for (int jj=0; jj<n_pmiss_cut; jj++){
      double fneg = 0;
      double fpos = 0;
      double fnegO = 0;
      double fposO = 0;
      double den = 0;
      double fnegSub = 0;
      double fposSub = 0;
      double truNeg = 0;
      double truPos = 0;
 
      for (int kk=0;kk<counter; kk++){
	/////////////////////////////////////////////
	//Calc using the smeared passing denominator
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//count smeared events that pass total (p/q, pmiss, mmiss, theta_pq)
	if (pmissSm[kk]>pmiss_cut[jj] && pmissSm[kk]<1.0 && mmissSm[kk]<mmiss_val[ii] && p_over_q_sm[kk] >= pq_min && p_over_q_sm[kk] <= pq_max_n){
	  den++;
	  //false positive: non-smeared fails, smeared passes
	  if (pmissO[kk]<pmiss_nomCut || mmissO[kk]>mmiss_nomCut || p_over_q_O[kk] < pq_min|| p_over_q_O[kk] > pq_max_n || pmissO[kk]>1){
	    fposO++;
	  }	  
	}
	//false negative: smeared fails, non-smeared passes
	if((pmissO[kk]>pmiss_nomCut && pmissO[kk]<1 && mmissO[kk]<mmiss_nomCut && p_over_q_O[kk] >= pq_min && p_over_q_O[kk] <= pq_max_n)
	   &&  (pmissSm[kk]<pmiss_cut[jj] || pmissSm[kk]>1 || mmissSm[kk]>mmiss_val[ii] || p_over_q_sm[kk] < pq_min || p_over_q_sm[kk] > pq_max_n)){
	  fnegO++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////
	//Calc using different denominators
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//false pos
	if ((pmissSm[kk]>pmiss_cut[jj] && pmissSm[kk]<1.0 && mmissSm[kk]<mmiss_val[ii] && p_over_q_sm[kk] >= pq_min && p_over_q_sm[kk] <= pq_max_n)
	    &&(pmissO[kk]<pmiss_nomCut || mmissO[kk]>mmiss_nomCut || p_over_q_O[kk] < pq_min|| p_over_q_O[kk] > pq_max_n || pmissO[kk]>1)){
	  fpos++;
	}

	//true negative
	if (pmissO[kk]<pmiss_nomCut || mmissO[kk]>mmiss_nomCut || p_over_q_O[kk] < pq_min|| p_over_q_O[kk] > pq_max_n || pmissO[kk]>1){
	  truNeg++;
	}

	//false neg
	if ((pmissO[kk]>pmiss_nomCut && pmissO[kk]<1 && mmissO[kk]<mmiss_nomCut && p_over_q_O[kk] >= pq_min && p_over_q_O[kk] <= pq_max_n)
	    &&  (pmissSm[kk]<pmiss_cut[jj] || pmissSm[kk]>1 || mmissSm[kk]>mmiss_val[ii] || p_over_q_sm[kk] < pq_min || p_over_q_sm[kk] > pq_max_n)){
	  fneg++;
	}

	//true positive
	if (pmissO[kk]>pmiss_nomCut && pmissO[kk]<1 && mmissO[kk]<mmiss_nomCut && p_over_q_O[kk] >= pq_min && p_over_q_O[kk] <= pq_max_n){
	  truPos++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
      }//end kk
      
      //false +ive: events that should not pass do pass
      falsePos[ii][jj] = fpos*100.0/(fpos+truNeg);
      falsePosO[ii][jj] = fposO*100.0/den;
      falsePosErr[ii][jj] = sqrt(fposO+pow(falsePosO[ii][jj]/100.0,2.0)*den)/den;
      //false -ive: events that should pass do not pass
      falseNeg[ii][jj] = fneg*100/(truPos+fneg);
      falseNegO[ii][jj] = fnegO*100/den;
      falseNegErr[ii][jj] = sqrt(fnegO+pow(falseNegO[ii][jj]/100.0,2.0)*den)/den;
      double diff = falseNegO[ii][jj] - falsePosO[ii][jj];
      double errDiff = sqrt(pow(falsePosErr[ii][jj],2.0)+pow(falseNegErr[ii][jj],2.0))*100.0;
      
      cout<<pmiss_cut[jj]<<","<<mmiss_val[ii]<<","<<falseNegO[ii][jj]<<","<<falsePosO[ii][jj]<<","<<falseNegErr[ii][jj]<<","<<falsePosErr[ii][jj]<<","<<diff<<","<<errDiff<<endl;

      if (abs(falsePosO[ii][jj]-falseNegO[ii][jj])<minimum){
	minimum = abs(falsePosO[ii][jj]-falseNegO[ii][jj]);
	mmissCut = ii;
	pmissCut = jj;
      }
      
    }//end jj
  }//end ii
  cout<<endl;
  //cout<<"optimal pmiss cut: "<<pmiss_cut[pmissCut]<<" optimal mmiss cut: "<<mmiss_val[mmissCut]<<" false Pos: "<<falsePosO[mmissCut][pmissCut]<<" false Neg: "<<falseNegO[mmissCut][pmissCut]<<endl;

  //false pos and false neg plots with different denominators
  TGraph *grP[n_mmiss_cut];
  TGraph *grN[n_mmiss_cut];
  TMultiGraph *mgP = new TMultiGraph();
  TMultiGraph *mgN = new TMultiGraph();
  TLegend* legendP = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendN = new TLegend(0.75,0.6,0.85,0.8);
  legendP->SetHeader("m_{miss}");
  legendP->SetBorderSize(0);
  legendN->SetHeader("m_{miss}");
  legendN->SetBorderSize(0);
  mgP->SetTitle("False Positive Rate;P_{miss} lower bound;False Positive [%]");
  mgN->SetTitle("False Negative Rate;P_{miss} lower bound;False Negative [%]");
  for (int ii=0; ii<n_mmiss_cut; ii++){
    grP[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falsePos[ii]);
    grP[ii]->SetMarkerStyle(20+ii);
    grP[ii]->SetMarkerColor(1+ii);
    mgP->Add(grP[ii],"p");
    legendP->AddEntry(grP[ii],Form("%.3f",mmiss_val[ii]),"p");
  
    grN[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falseNeg[ii]);
    grN[ii]->SetMarkerStyle(20+ii);
    grN[ii]->SetMarkerColor(1+ii);
    mgN->Add(grN[ii],"p");
    legendN->AddEntry(grN[ii],Form("%.3f",mmiss_val[ii]),"p");
  }

  //false pos and false neg plot with smeared passing denominator
  TGraph *grPO[n_mmiss_cut];
  TGraph *grNO[n_mmiss_cut];
  TMultiGraph *mgPO = new TMultiGraph();
  TMultiGraph *mgNO = new TMultiGraph();
  TLegend* legendPO = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendNO = new TLegend(0.75,0.6,0.85,0.8);
  legendPO->SetHeader("m_{miss}");
  legendPO->SetBorderSize(0);
  legendNO->SetHeader("m_{miss}");
  legendNO->SetBorderSize(0);
  mgPO->SetTitle("False Positive Rate;P_{miss} lower bound;False Positive Rate [%]");
  mgNO->SetTitle("False Negative per smeared passing;P_{miss} lower bound;False Negative per smeared passing [%]");
  for (int ii=0; ii<n_mmiss_cut; ii++){
    grPO[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falsePosO[ii]);
    grPO[ii]->SetMarkerStyle(20+ii);
    grPO[ii]->SetMarkerColor(1+ii);
    mgPO->Add(grPO[ii],"p");
    legendPO->AddEntry(grPO[ii],Form("%.3f",mmiss_val[ii]),"p");
  
    grNO[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falseNegO[ii]);
    grNO[ii]->SetMarkerStyle(20+ii);
    grNO[ii]->SetMarkerColor(1+ii);
    mgNO->Add(grNO[ii],"p");
    legendNO->AddEntry(grNO[ii],Form("%.3f",mmiss_val[ii]),"p");
  }

  //f->Close();
  double totEvents;
  canvas->Update();

  mgP->Draw("a");
  mgP->GetYaxis()->SetRangeUser(0,100);
  legendP->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgN->Draw("a");
  mgN->GetYaxis()->SetRangeUser(0,100);
  legendN->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgPO->Draw("a");
  mgPO->GetYaxis()->SetRangeUser(0,100);
  legendPO->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgNO->Draw("a");
  mgNO->GetYaxis()->SetRangeUser(0,100);
  legendNO->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
  hfid_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  hfid_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  hec_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  hec_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_etheta->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emom->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pmiss->SetLineColor(kRed);
  h_pmissSmear->SetLineColor(kBlue);
  h_pmiss->Draw();
  h_pmissSmear->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pmissCut->SetLineColor(kBlue);
  eg2_pmiss->SetLineColor(kMagenta);
  eg2_pmiss->Scale(h_pmissCut->GetMaximum()/eg2_pmiss->GetMaximum());
  h_pmissCut->Draw();
  eg2_pmiss->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emissCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  
  h_emiss->SetLineColor(kRed);
  h_emissSmear->SetLineColor(kBlue);
  h_emiss->Draw();
  h_emissSmear->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmomCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom->SetLineColor(kRed);
  h_nmomSmear->SetLineColor(kBlue);
  h_nmom->Draw();
  h_nmomSmear->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nthetaCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_ntheta->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapqCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_q2->Draw();
  totEvents = h_q2->GetEntries();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_smearCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_smear->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hist_pq_thetapq->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  hist_pq_thetapq_smear->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  hist_pmiss_xb->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpmSm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmissCut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_mmiss->SetLineColor(kRed);
  h_mmissSmear->SetLineColor(kBlue);
  h_mmiss->Draw();
  h_mmissSmear->Draw("same");
  //h_zvtx_precut->Scale(1./h_zvtx_precut->GetMaximum());
  //h_zvtx_cut->Scale(1./h_zvtx_cut->GetMaximum());
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_zvtx_precut->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_timeres->Draw();
  TF1 *g1 = new TF1("g1","gaus",-1,1);
  h_timeres->Fit(g1,"R");
  g1->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pmMpm->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  pvps->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmvpm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmvpmS->Draw("colz");
  canvas->Print( (pdf_file_name + ")").c_str());



  fout->Write();
  fout->Close();
  cout<<"Good events: "<<totEvents<<endl;

}
