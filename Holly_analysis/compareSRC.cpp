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
const double BE_2H   = 0.00222;			// 2H  binding energy in GeV
const double BE_3H   = 0.00848;			// 3H  binding energy in GeV
const double BE_3He  = 0.00772;                 // 3He binding energy in GeV
const double BE_4He  = 0.02830;                 // 4He binding energy in GeV
const double BE_11B  = 0.07620;			// 11B binding energy in GeV
const double BE_11C  = 0.07344;			// 11C binding energy in GeV
const double BE_12C  = 0.09215;                 // 12C binding energy in GeV
const double M2H     =   Mp +   Mn -BE_2H ;
const double M3H     =   Mp + 2*Mn -BE_3H ;
const double M3He    = 2* Mp+Mn -BE_3He;        // 3He mass in GeV
const double M4He    = 2*(Mp+Mn)-BE_4He;        // 4He mass in GeV
const double M11B    = 5*Mp + 6*Mn -BE_11B;
const double M11C    = 6*Mp + 5*Mn -BE_11C;
const double M12C    = 6*(Mp+Mn)-BE_12C;        // 12C mass in GeV
const double pi      = 3.14159;
const double cc      = 2.99792458E+8;

//vertex cuts
const double vtx_min_C = 4.0;//cm
const double vtx_max_C = 7.0;//cm
const double vtx_min_He = -2.5;//cm
const double vtx_max_He = -0.5;//cm

//SRC cuts defined:
const double xbP_cut = 1.2;
const double xbN_cut = 1.1;
const double thetapqSRC_cut = 25.0;//deg
const double pq_min = 0.62;
const double pq_max_n = 1.1;
const double pq_max_p = 0.96;
const double mmiss_nomCut_n = 1.175;
const double pmissSRC_nomCut_n = 0.38;
const double mmiss_nomCut_p = 1.1;
const double pmissSRC_nomCut_p = 0.3;

//MF cuts defined:
const double y_cut_min = -0.05;
const double y_cut_max = 0.2;
const double nu_min = 0.9;
const double nu_max = 1.6;
const double thetapqMF_cut = 7.0;//deg
const double emiss_nomCut = 0.16;
const double pmissMF_nomCut = 0.325;


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

        double ec_extr_res = -0.09;//0967682;
	double dt = 0.392319;//0.625                 //ec_time_res;
	const double d      = 516.9152;       //cm, distance from interaction vtx to EC
	const double Mn     = 0.939565378;    // Neutron mass in GeV
	const double c_cm_s = 2.99792458E+10; // cm/s, speed of light
	const double ns_2_s = 1E-09;

	double val = p*p/(Mn*Mn)*sqrt(Mn*Mn+p*p)*dt/d*c_cm_s*ns_2_s;
	return sqrt(val*val/p/p+ec_extr_res*ec_extr_res)*p;
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

int main(){

  const int maxParticles=50;
  TFile *fout = new TFile("output/out_SRC.root","RECREATE");

  TCanvas *canvas = new TCanvas("canvas","output plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/out_SRC.pdf";
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  //////////////////////////////////
  // Plots/////////////////////////
  /////////////////////////////////
  TH1F* h_emom[3][3]; 
  TH1F* h_etheta[3][3];
  TH1F* h_q2[3][3];
  TH1F* h_nu[3][3];
  TH1F* h_xb[3][3];
  TH1F* h_pmiss[3][3];
  TH2F* h_emVpm[3][3];
  TH1F* h_thetapq[3][3];
  TH1F* h_nmom[3][3];
  TH1F* h_ntheta[3][3];
  TH1D * h_emiss[3][3]; 
  TH1D * h_mmiss[3][3];


  //target: 3He, 4He, 12C
  for (int ii=0;ii<3;ii++){
    //lead: proton, smeared proton, neutron
    for(int jj=0;jj<3;jj++){
      h_emom[ii][jj]=new TH1F(Form("h_emom_%i_%i",ii,jj),";P_{e-} [GeV/c] (with lead cuts);",50,0,6); 
      h_etheta[ii][jj]=new TH1F(Form("h_etheta_%i_%i",ii,jj),";#theta_{e-} [deg] (with lead cuts);",50,0,45);
      h_q2[ii][jj]=new TH1F(Form("h_q2_%i_%i",ii,jj),";Q^{2} (with lead cuts);",50,0,6);
      h_nu[ii][jj]=new TH1F(Form("h_nu_%i_%i",ii,jj),";#nu [GeV] (with lead cuts);",50,0,3);
      h_xb[ii][jj]=new TH1F(Form("h_xb_%i_%i",ii,jj),";x_{B} (with lead cuts);",50,0,2);
      h_pmiss[ii][jj]=new TH1F(Form("h_pmiss_%i_%i",ii,jj),";P_{miss} [GeV/c];",50,0,1.5);
      h_emVpm[ii][jj]=new TH2F(Form("h_emVpm_%i_%i",ii,jj),";pmiss;emiss",50,0,0.8,50,-0.5,0.5);
      h_thetapq[ii][jj]=new TH1F(Form("h_thetapq_%i_%i",ii,jj),";#theta_{lead N/q} [deg];",50,0,25);
      h_nmom[ii][jj]=new TH1F(Form("h_nmom_%i_%i",ii,jj),";P_{lead N} [GeV/c];",50,0,3);
      h_ntheta[ii][jj]=new TH1F(Form("h_ntheta_%i_%i",ii,jj),";#theta_{lead N} [deg];",50,0,90);
      h_emiss[ii][jj]= new TH1D(Form("h_emiss_%i_%i",ii,jj),";missing energy",50,-1,2);
      h_mmiss[ii][jj]= new TH1D(Form("h_mmiss_%i_%i",ii,jj),";missing mass [GeV/c^{2}]",50,0,3);
    }
  }
  
 
  //////////////////////////////////
  // Set some target dep parameters
  /////////////////////////////////
  double tgt_min, tgt_max, tgt_M, tgt_R;
  int tgt_Z, tgt_A;

  ////////////////////////////////
  //Helium-3
  ////////////////////////////////
  TString target = "3He";
  tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 3; tgt_M = M3He; tgt_R = M2H;
  Fiducial fid_params(4461,2250,5996,"3He",true);

  //open the 3He proton file
  TFile *f_3p = new TFile("../../../data/19jun20/skim_He3_p_Jun20.root");
  TTree * intree = (TTree*)f_3p->Get("T");
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles];
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  gRandom = new TRandom3();
  double theta_prev = 0;
  double p_over_q_prev = 0;

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
      if (xB<xbN_cut){continue;}

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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  	  
	  //smear the proton momentum resolution
	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
	  TVector3 u1 = mom_lead.Unit();
	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

	  //recalculate the smeared factors
	  TVector3 mom_missSmear = mom_smear - q;
	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  if(mom_smear.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();
	  double p_over_q_smear = mom_smear.Mag() / q.Mag();



	  //select the 3He proton sample
	  if (p_over_q>pq_min && p_over_q<pq_max_p && mmiss<mmiss_nomCut_p && xB>xbP_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_p){
	    h_emom[0][0]->Fill(mom_e.Mag()); 
	    h_etheta[0][0]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[0][0]->Fill(QSq);
	    h_nu[0][0]->Fill(Nu);
	    h_xb[0][0]->Fill(xB);
	    h_pmiss[0][0]->Fill(mom_miss.Mag());
	    h_emVpm[0][0]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[0][0]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[0][0]->Fill(mom_lead.Mag());
	    h_ntheta[0][0]->Fill(theta);
	    h_emiss[0][0]->Fill(emiss);
	    h_mmiss[0][0]->Fill(mmiss);	    
	  }

	  //select the 3He smeared proton sample
	  if (p_over_q_smear>pq_min && p_over_q_smear<pq_max_n && mmissSmear<mmiss_nomCut_n && xB>xbN_cut && mom_missSmear.Mag()<1 && mom_missSmear.Mag()>pmissSRC_nomCut_n){
	    h_emom[0][1]->Fill(mom_e.Mag()); 
	    h_etheta[0][1]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[0][1]->Fill(QSq);
	    h_nu[0][1]->Fill(Nu);
	    h_xb[0][1]->Fill(xB);
	    h_pmiss[0][1]->Fill(mom_missSmear.Mag());
	    h_emVpm[0][1]->Fill(mom_missSmear.Mag(),emissSmear);
	    h_thetapq[0][1]->Fill(mom_smear.Angle(q)*180./M_PI);
	    h_nmom[0][1]->Fill(mom_smear.Mag());
	    h_ntheta[0][1]->Fill(theta);
	    h_emiss[0][1]->Fill(emissSmear);
	    h_mmiss[0][1]->Fill(mmissSmear);
	    
	  }
	  

	}//end leading
    }//end event
  
  intree->Delete();
  
  //open the 3He neutron file
  TFile *f_3n = new TFile("../../../data/19jun20/skim_He3_n_Jun20.root");
  intree = (TTree*)f_3n->Get("T");
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  int nentries = intree->GetEntries();
  cout<<"number of n entries: "<<nentries<<endl;
  
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
      if (xB<xbN_cut){continue;}

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at neutrons
	  if (type[part] != 2112)
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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading neutron!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();


	  //select the 3He neutron sample
	  if (p_over_q>pq_min && p_over_q<pq_max_n && mmiss<mmiss_nomCut_n && xB>xbN_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_n){

	    h_emom[0][2]->Fill(mom_e.Mag()); 
	    h_etheta[0][2]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[0][2]->Fill(QSq);
	    h_nu[0][2]->Fill(Nu);
	    h_xb[0][2]->Fill(xB);
	    h_pmiss[0][2]->Fill(mom_miss.Mag());
	    h_emVpm[0][2]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[0][2]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[0][2]->Fill(mom_lead.Mag());
	    h_ntheta[0][2]->Fill(theta);
	    h_emiss[0][2]->Fill(emiss);
	    h_mmiss[0][2]->Fill(mmiss);
	    
	  }	  

	}//end leading
    }//end event
    
  
  ////////////////////////////////
  //Helium-4
  ////////////////////////////////
  intree->Delete();
  
  target = "4He";
  tgt_min = tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 4; tgt_M = M4He; tgt_R = M3H;
  Fiducial fid_params4He(4461,2250,5996,"4He",true);
  //open the 4He proton file
  TFile *f_4p = new TFile("../../../data/19jun20/skim_He4_p_Jun20.root");
  intree = (TTree*)f_4p->Get("T");
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

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
      if (xB<xbN_cut){continue;}

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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  	  
	  //smear the proton momentum resolution
	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
	  TVector3 u1 = mom_lead.Unit();
	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

	  //recalculate the smeared factors
	  TVector3 mom_missSmear = mom_smear - q;
	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params4He.pFiducialCut(mom_lead)) {continue;}
	  
	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  if(mom_smear.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();
	  double p_over_q_smear = mom_smear.Mag() / q.Mag();



	  //select the 4He proton sample
	  if (p_over_q>pq_min && p_over_q<pq_max_p && mmiss<mmiss_nomCut_p && xB>xbP_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_p){
	    h_emom[1][0]->Fill(mom_e.Mag()); 
	    h_etheta[1][0]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[1][0]->Fill(QSq);
	    h_nu[1][0]->Fill(Nu);
	    h_xb[1][0]->Fill(xB);
	    h_pmiss[1][0]->Fill(mom_miss.Mag());
	    h_emVpm[1][0]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[1][0]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[1][0]->Fill(mom_lead.Mag());
	    h_ntheta[1][0]->Fill(theta);
	    h_emiss[1][0]->Fill(emiss);
	    h_mmiss[1][0]->Fill(mmiss);	    
	  }

	  //select the 4He smeared proton sample
	  if (p_over_q_smear>pq_min && p_over_q_smear<pq_max_n && mmissSmear<mmiss_nomCut_n && xB>xbN_cut && mom_missSmear.Mag()<1 && mom_missSmear.Mag()>pmissSRC_nomCut_n){
	    h_emom[1][1]->Fill(mom_e.Mag()); 
	    h_etheta[1][1]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[1][1]->Fill(QSq);
	    h_nu[1][1]->Fill(Nu);
	    h_xb[1][1]->Fill(xB);
	    h_pmiss[1][1]->Fill(mom_missSmear.Mag());
	    h_emVpm[1][1]->Fill(mom_missSmear.Mag(),emissSmear);
	    h_thetapq[1][1]->Fill(mom_smear.Angle(q)*180./M_PI);
	    h_nmom[1][1]->Fill(mom_smear.Mag());
	    h_ntheta[1][1]->Fill(theta);
	    h_emiss[1][1]->Fill(emissSmear);
	    h_mmiss[1][1]->Fill(mmissSmear);
	    
	  }
	  

	}//end leading
    }//end event
  
  intree->Delete();
  //open the 4He neutron file
  TFile *f_4n = new TFile("../../../data/19jun20/skim_He4_n_Jun20.root");
  intree = (TTree*)f_4n->Get("T");
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
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
      if (xB<xbN_cut){continue;}

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at neutrons
	  if (type[part] != 2112)
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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading neutron!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params4He.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();


	  //select the 4He neutron sample
	  if (p_over_q>pq_min && p_over_q<pq_max_n && mmiss<mmiss_nomCut_n && xB>xbN_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_n){

	    h_emom[1][2]->Fill(mom_e.Mag()); 
	    h_etheta[1][2]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[1][2]->Fill(QSq);
	    h_nu[1][2]->Fill(Nu);
	    h_xb[1][2]->Fill(xB);
	    h_pmiss[1][2]->Fill(mom_miss.Mag());
	    h_emVpm[1][2]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[1][2]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[1][2]->Fill(mom_lead.Mag());
	    h_ntheta[1][2]->Fill(theta);
	    h_emiss[1][2]->Fill(emiss);
	    h_mmiss[1][2]->Fill(mmiss);
	    
	  }	  

	}//end leading
    }//end event
    

  
  ////////////////////////////////
  //Carbon-12
  ////////////////////////////////
  target = "12C";
  tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C; tgt_R = M11B;
  Fiducial fid_paramsC(4461,2250,5996,"12C",true);
  intree->Delete();
  
  //open the 12C proton file
  TFile *f_5p = new TFile("../../../data/19jun20/skim_C12_p_Jun20.root");
  intree = (TTree*)f_5p->Get("T");
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

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
      if (xB<xbN_cut){continue;}

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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  	  
	  //smear the proton momentum resolution
	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
	  TVector3 u1 = mom_lead.Unit();
	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

	  //recalculate the smeared factors
	  TVector3 mom_missSmear = mom_smear - q;
	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_paramsC.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  if(mom_smear.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();
	  double p_over_q_smear = mom_smear.Mag() / q.Mag();



	  //select the 12C proton sample
	  if (p_over_q>pq_min && p_over_q<pq_max_p && mmiss<mmiss_nomCut_p && xB>xbP_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_p){
	    h_emom[2][0]->Fill(mom_e.Mag()); 
	    h_etheta[2][0]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[2][0]->Fill(QSq);
	    h_nu[2][0]->Fill(Nu);
	    h_xb[2][0]->Fill(xB);
	    h_pmiss[2][0]->Fill(mom_miss.Mag());
	    h_emVpm[2][0]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[2][0]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[2][0]->Fill(mom_lead.Mag());
	    h_ntheta[2][0]->Fill(theta);
	    h_emiss[2][0]->Fill(emiss);
	    h_mmiss[2][0]->Fill(mmiss);	    
	  }

	  //select the 12C smeared proton sample
	  if (p_over_q_smear>pq_min && p_over_q_smear<pq_max_n && mmissSmear<mmiss_nomCut_n && xB>xbN_cut && mom_missSmear.Mag()<1 && mom_missSmear.Mag()>pmissSRC_nomCut_n){
	    h_emom[2][1]->Fill(mom_e.Mag()); 
	    h_etheta[2][1]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[2][1]->Fill(QSq);
	    h_nu[2][1]->Fill(Nu);
	    h_xb[2][1]->Fill(xB);
	    h_pmiss[2][1]->Fill(mom_missSmear.Mag());
	    h_emVpm[2][1]->Fill(mom_missSmear.Mag(),emissSmear);
	    h_thetapq[2][1]->Fill(mom_smear.Angle(q)*180./M_PI);
	    h_nmom[2][1]->Fill(mom_smear.Mag());
	    h_ntheta[2][1]->Fill(theta);
	    h_emiss[2][1]->Fill(emissSmear);
	    h_mmiss[2][1]->Fill(mmissSmear);
	    
	  }
	  

	}//end leading
    }//end event
  
  intree->Delete();
  //open the 12C neutron file
  TFile *f_5n = new TFile("../../../data/19jun20/skim_C12_n_Jun20.root");
  intree = (TTree*)f_5n->Get("T");
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
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
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
      if (xB<xbN_cut){continue;}

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at neutrons
	  if (type[part] != 2112)
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
	  if (theta_pq> thetapqSRC_cut){continue;}

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
	  else if (abs(p_over_q-pq_min)/pq_min*abs(theta_pq-thetapqSRC_cut)/thetapqSRC_cut < abs(p_over_q_prev-pq_min)/pq_min*abs(theta_prev-thetapqSRC_cut)/thetapqSRC_cut)
	    {leadingID=part;}
	  	  
	}

      //Now we have a leading neutron!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_paramsC.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;

	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();


	  //select the 12C neutron sample
	  if (p_over_q>pq_min && p_over_q<pq_max_n && mmiss<mmiss_nomCut_n && xB>xbN_cut && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_n){

	    h_emom[2][2]->Fill(mom_e.Mag()); 
	    h_etheta[2][2]->Fill(mom_e.Theta()*180./M_PI);
	    h_q2[2][2]->Fill(QSq);
	    h_nu[2][2]->Fill(Nu);
	    h_xb[2][2]->Fill(xB);
	    h_pmiss[2][2]->Fill(mom_miss.Mag());
	    h_emVpm[2][2]->Fill(mom_miss.Mag(),emiss);
	    h_thetapq[2][2]->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_nmom[2][2]->Fill(mom_lead.Mag());
	    h_ntheta[2][2]->Fill(theta);
	    h_emiss[2][2]->Fill(emiss);
	    h_mmiss[2][2]->Fill(mmiss);
	    
	  }	  

	}//end leading
    }//end event
  
  //f->Close();
  int protons_3He, protons_4He, protons_12C, smeared_3He, smeared_4He, smeared_12C, neutrons_3He, neutrons_4He, neutrons_12C;
  protons_3He = h_emom[0][0]->GetEntries();
  neutrons_3He = h_emom[0][2]->GetEntries();
  smeared_3He =  h_emom[0][1]->GetEntries();
  protons_4He = h_emom[1][0]->GetEntries();
  neutrons_4He = h_emom[1][2]->GetEntries();
  smeared_4He =  h_emom[1][1]->GetEntries();
  protons_12C = h_emom[2][0]->GetEntries();
  neutrons_12C = h_emom[2][2]->GetEntries();
  smeared_12C =  h_emom[2][1]->GetEntries();
  
  canvas->Update();

  h_emVpm[0][0]->SetTitle("3He protons");
  h_emVpm[0][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[0][1]->SetTitle("3He smeared protons");
  h_emVpm[0][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[0][2]->SetTitle("3He neutrons");
  h_emVpm[0][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm[1][0]->SetTitle("4He protons");
  h_emVpm[1][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[1][1]->SetTitle("4He smeared protons");
  h_emVpm[1][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[1][2]->SetTitle("4He neutrons");
  h_emVpm[1][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm[2][0]->SetTitle("12C protons");
  h_emVpm[2][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[2][1]->SetTitle("12C smeared protons");
  h_emVpm[2][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[2][2]->SetTitle("12C neutrons");
  h_emVpm[2][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  ///////////////////////////////////////
  //Compare neutrons and smeared protons
  ///////////////////////////////////////
  //Helium-3
  h_emom[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_emom[0][1]->SetLineColor(kRed);
  h_emom[0][2]->SetLineColor(kBlue);
  h_emom[0][1]->Scale(1./h_emom[0][1]->Integral());
  h_emom[0][2]->Scale(1./h_emom[0][2]->Integral());
  h_emom[0][1]->Draw("hist");
  h_emom[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_etheta[0][1]->SetLineColor(kRed);
  h_etheta[0][2]->SetLineColor(kBlue);
  h_etheta[0][1]->Scale(1./h_etheta[0][1]->Integral());
  h_etheta[0][2]->Scale(1./h_etheta[0][2]->Integral());
  h_etheta[0][1]->Draw("hist");
  h_etheta[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_q2[0][1]->SetLineColor(kRed);
  h_q2[0][2]->SetLineColor(kBlue);
  h_q2[0][1]->Scale(1./h_q2[0][1]->Integral());
  h_q2[0][2]->Scale(1./h_q2[0][2]->Integral());
  h_q2[0][1]->Draw("hist");
  h_q2[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_nu[0][1]->SetLineColor(kRed);
  h_nu[0][2]->SetLineColor(kBlue);
  h_nu[0][1]->Scale(1./h_nu[0][1]->Integral());
  h_nu[0][2]->Scale(1./h_nu[0][2]->Integral());
  h_nu[0][1]->Draw("hist");
  h_nu[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_xb[0][1]->SetLineColor(kRed);
  h_xb[0][2]->SetLineColor(kBlue);
  h_xb[0][1]->Scale(1./h_xb[0][1]->Integral());
  h_xb[0][2]->Scale(1./h_xb[0][2]->Integral());
  h_xb[0][1]->Draw("hist");
  h_xb[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_pmiss[0][1]->SetLineColor(kRed);
  h_pmiss[0][2]->SetLineColor(kBlue);
  h_pmiss[0][1]->Scale(1./h_pmiss[0][1]->Integral());
  h_pmiss[0][2]->Scale(1./h_pmiss[0][2]->Integral());
  h_pmiss[0][1]->Draw("hist");
  h_pmiss[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_thetapq[0][1]->SetLineColor(kRed);
  h_thetapq[0][2]->SetLineColor(kBlue);
  h_thetapq[0][1]->Scale(1./h_thetapq[0][1]->Integral());
  h_thetapq[0][2]->Scale(1./h_thetapq[0][2]->Integral());
  h_thetapq[0][1]->Draw("hist");
  h_thetapq[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_nmom[0][1]->SetLineColor(kRed);
  h_nmom[0][2]->SetLineColor(kBlue);
  h_nmom[0][1]->Scale(1./h_nmom[0][1]->Integral());
  h_nmom[0][2]->Scale(1./h_nmom[0][2]->Integral());
  h_nmom[0][1]->Draw("hist");
  h_nmom[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_ntheta[0][1]->SetLineColor(kRed);
  h_ntheta[0][2]->SetLineColor(kBlue);
  h_ntheta[0][1]->Scale(1./h_ntheta[0][1]->Integral());
  h_ntheta[0][2]->Scale(1./h_ntheta[0][2]->Integral());
  h_ntheta[0][1]->Draw("hist");
  h_ntheta[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_emiss[0][1]->SetLineColor(kRed);
  h_emiss[0][2]->SetLineColor(kBlue);
  h_emiss[0][1]->Scale(1./h_emiss[0][1]->Integral());
  h_emiss[0][2]->Scale(1./h_emiss[0][2]->Integral());
  h_emiss[0][1]->Draw("hist");
  h_emiss[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_mmiss[0][1]->SetLineColor(kRed);
  h_mmiss[0][2]->SetLineColor(kBlue);
  h_mmiss[0][1]->Scale(1./h_mmiss[0][1]->Integral());
  h_mmiss[0][2]->Scale(1./h_mmiss[0][2]->Integral());
  h_mmiss[0][1]->Draw("hist");
  h_mmiss[0][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  //Helium-4
  h_emom[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_emom[1][1]->SetLineColor(kRed);
  h_emom[1][2]->SetLineColor(kBlue);
  h_emom[1][1]->Scale(1./h_emom[1][1]->Integral());
  h_emom[1][2]->Scale(1./h_emom[1][2]->Integral());
  h_emom[1][1]->Draw("hist");
  h_emom[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_etheta[1][1]->SetLineColor(kRed);
  h_etheta[1][2]->SetLineColor(kBlue);
  h_etheta[1][1]->Scale(1./h_etheta[1][1]->Integral());
  h_etheta[1][2]->Scale(1./h_etheta[1][2]->Integral());
  h_etheta[1][1]->Draw("hist");
  h_etheta[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_q2[1][1]->SetLineColor(kRed);
  h_q2[1][2]->SetLineColor(kBlue);
  h_q2[1][1]->Scale(1./h_q2[1][1]->Integral());
  h_q2[1][2]->Scale(1./h_q2[1][2]->Integral());
  h_q2[1][1]->Draw("hist");
  h_q2[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_nu[1][1]->SetLineColor(kRed);
  h_nu[1][2]->SetLineColor(kBlue);
  h_nu[1][1]->Scale(1./h_nu[1][1]->Integral());
  h_nu[1][2]->Scale(1./h_nu[1][2]->Integral());
  h_nu[1][1]->Draw("hist");
  h_nu[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_xb[1][1]->SetLineColor(kRed);
  h_xb[1][2]->SetLineColor(kBlue);
  h_xb[1][1]->Scale(1./h_xb[1][1]->Integral());
  h_xb[1][2]->Scale(1./h_xb[1][2]->Integral());
  h_xb[1][1]->Draw("hist");
  h_xb[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_pmiss[1][1]->SetLineColor(kRed);
  h_pmiss[1][2]->SetLineColor(kBlue);
  h_pmiss[1][1]->Scale(1./h_pmiss[1][1]->Integral());
  h_pmiss[1][2]->Scale(1./h_pmiss[1][2]->Integral());
  h_pmiss[1][1]->Draw("hist");
  h_pmiss[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_thetapq[1][1]->SetLineColor(kRed);
  h_thetapq[1][2]->SetLineColor(kBlue);
  h_thetapq[1][1]->Scale(1./h_thetapq[1][1]->Integral());
  h_thetapq[1][2]->Scale(1./h_thetapq[1][2]->Integral());
  h_thetapq[1][1]->Draw("hist");
  h_thetapq[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_nmom[1][1]->SetLineColor(kRed);
  h_nmom[1][2]->SetLineColor(kBlue);
  h_nmom[1][1]->Scale(1./h_nmom[1][1]->Integral());
  h_nmom[1][2]->Scale(1./h_nmom[1][2]->Integral());
  h_nmom[1][1]->Draw("hist");
  h_nmom[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_ntheta[1][1]->SetLineColor(kRed);
  h_ntheta[1][2]->SetLineColor(kBlue);
  h_ntheta[1][1]->Scale(1./h_ntheta[1][1]->Integral());
  h_ntheta[1][2]->Scale(1./h_ntheta[1][2]->Integral());
  h_ntheta[1][1]->Draw("hist");
  h_ntheta[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_emiss[1][1]->SetLineColor(kRed);
  h_emiss[1][2]->SetLineColor(kBlue);
  h_emiss[1][1]->Scale(1./h_emiss[1][1]->Integral());
  h_emiss[1][2]->Scale(1./h_emiss[1][2]->Integral());
  h_emiss[1][1]->Draw("hist");
  h_emiss[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_mmiss[1][1]->SetLineColor(kRed);
  h_mmiss[1][2]->SetLineColor(kBlue);
  h_mmiss[1][1]->Scale(1./h_mmiss[1][1]->Integral());
  h_mmiss[1][2]->Scale(1./h_mmiss[1][2]->Integral());
  h_mmiss[1][1]->Draw("hist");
  h_mmiss[1][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  //Carbon-12
  h_emom[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_emom[2][1]->SetLineColor(kRed);
  h_emom[2][2]->SetLineColor(kBlue);
  h_emom[2][1]->Scale(1./h_emom[2][1]->Integral());
  h_emom[2][2]->Scale(1./h_emom[2][2]->Integral());
  h_emom[2][1]->Draw("hist");
  h_emom[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_etheta[2][1]->SetLineColor(kRed);
  h_etheta[2][2]->SetLineColor(kBlue);
  h_etheta[2][1]->Scale(1./h_etheta[2][1]->Integral());
  h_etheta[2][2]->Scale(1./h_etheta[2][2]->Integral());
  h_etheta[2][1]->Draw("hist");
  h_etheta[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_q2[2][1]->SetLineColor(kRed);
  h_q2[2][2]->SetLineColor(kBlue);
  h_q2[2][1]->Scale(1./h_q2[2][1]->Integral());
  h_q2[2][2]->Scale(1./h_q2[2][2]->Integral());
  h_q2[2][1]->Draw("hist");
  h_q2[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_nu[2][1]->SetLineColor(kRed);
  h_nu[2][2]->SetLineColor(kBlue);
  h_nu[2][1]->Scale(1./h_nu[2][1]->Integral());
  h_nu[2][2]->Scale(1./h_nu[2][2]->Integral());
  h_nu[2][1]->Draw("hist");
  h_nu[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_xb[2][1]->SetLineColor(kRed);
  h_xb[2][2]->SetLineColor(kBlue);
  h_xb[2][1]->Scale(1./h_xb[2][1]->Integral());
  h_xb[2][2]->Scale(1./h_xb[2][2]->Integral());
  h_xb[2][1]->Draw("hist");
  h_xb[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_pmiss[2][1]->SetLineColor(kRed);
  h_pmiss[2][2]->SetLineColor(kBlue);
  h_pmiss[2][1]->Scale(1./h_pmiss[2][1]->Integral());
  h_pmiss[2][2]->Scale(1./h_pmiss[2][2]->Integral());
  h_pmiss[2][1]->Draw("hist");
  h_pmiss[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_thetapq[2][1]->SetLineColor(kRed);
  h_thetapq[2][2]->SetLineColor(kBlue);
  h_thetapq[2][1]->Scale(1./h_thetapq[2][1]->Integral());
  h_thetapq[2][2]->Scale(1./h_thetapq[2][2]->Integral());
  h_thetapq[2][1]->Draw("hist");
  h_thetapq[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_nmom[2][1]->SetLineColor(kRed);
  h_nmom[2][2]->SetLineColor(kBlue);
  h_nmom[2][1]->Scale(1./h_nmom[2][1]->Integral());
  h_nmom[2][2]->Scale(1./h_nmom[2][2]->Integral());
  h_nmom[2][1]->Draw("hist");
  h_nmom[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_ntheta[2][1]->SetLineColor(kRed);
  h_ntheta[2][2]->SetLineColor(kBlue);
  h_ntheta[2][1]->Scale(1./h_ntheta[2][1]->Integral());
  h_ntheta[2][2]->Scale(1./h_ntheta[2][2]->Integral());
  h_ntheta[2][1]->Draw("hist");
  h_ntheta[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_emiss[2][1]->SetLineColor(kRed);
  h_emiss[2][2]->SetLineColor(kBlue);
  h_emiss[2][1]->Scale(1./h_emiss[2][1]->Integral());
  h_emiss[2][2]->Scale(1./h_emiss[2][2]->Integral());
  h_emiss[2][1]->Draw("hist");
  h_emiss[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_mmiss[2][1]->SetLineColor(kRed);
  h_mmiss[2][2]->SetLineColor(kBlue);
  h_mmiss[2][1]->Scale(1./h_mmiss[2][1]->Integral());
  h_mmiss[2][2]->Scale(1./h_mmiss[2][2]->Integral());
  h_mmiss[2][1]->Draw("hist");
  h_mmiss[2][2]->Draw("histsame");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  //////////////////////////
  //Compare amongst targets
  //////////////////////////
  /*

  //h_emom[0][0]->SetTitle("protons (3He-r, 4He-b, 12C-p)");
  //h_emom[0][0]->SetLineColor(kRed);
  //h_emom[1][0]->SetLineColor(kBlue);
  //h_emom[2][0]->SetLineColor("kMagenta");
  //h_emom[0][0]->Scale(1./h_emom[0][0]->Integral());
  //h_emom[1][0]->Scale(1./h_emom[1][0]->Integral());
  //h_emom[2][0]->Scale(1./h_emom[2][0]->Integral());
  //h_emom[0][0]->Draw("hist");
  //h_emom[1][0]->Scale(1./h_emom[1][0]->Integral());
  //h_emom[2][0]->Scale(1./h_emom[2][0]->Integral());
  */
  canvas->Print( (pdf_file_name + ")").c_str());

  cout<<"3He protons: "<<protons_3He<<" 3He neutrons: "<<neutrons_3He<<" 3He smeared: "<<smeared_3He<<endl;
  cout<<"4He protons: "<<protons_4He<<" 4He neutrons: "<<neutrons_4He<<" 4He smeared: "<<smeared_4He<<endl;
  cout<<"12C protons: "<<protons_12C<<" 12C neutrons: "<<neutrons_12C<<" 12C smeared: "<<smeared_12C<<endl;

  

  fout->Write();
  fout->Close();

}
