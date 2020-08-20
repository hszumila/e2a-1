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

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  //plot theta vs phi before
  TH2F *hfid_pre = new TH2F("hfid_pre","no p fid cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot theta vs phi after
  TH2F *hfid_post = new TH2F("hfid_post","with p fid cut & EC 10cm cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot EC x vs y before
  TH2F *hec_pre = new TH2F("hec_pre","no fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
  //plot EC x vs y after
  TH2F *hec_post = new TH2F("hec_post","10cm fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
  //plot pmiss and xb
  TH2D * hist_pmiss_xb = new TH2D("pmiss_xb","Events before proton cuts;pmiss;Xb;Counts",28,0.3,1.0,40,1,3);
  //plot theta_pq and p/q
  TH2D * hist_pq_thetapq = new TH2D("hist_pq_thetapq","Events before proton cuts;p/q;#theta_{pq} [deg];Counts",100,0,1.5,100,0,70);
  //plot theta_pq and p/q
  TH2D * hist_pq_thetapq_smear = new TH2D("hist_pq_thetapq_smear","Smeared events before lead cuts;p_{smeared}/q;#theta_{pq} [deg];Counts",100,0,1.5,100,0,70);
  //plot z vertex
  TH1D * h_zvtx_precut = new TH1D("h_zvtx_precut","Corrected z vertex-no cut;z vertex [cm]",100,-10,10);
  //plot the missing mass
  TH1D * h_mmiss = new TH1D("h_mmiss","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  TH1D * h_mmissSmear = new TH1D("h_mmissSmear","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  TH1D * h_mmissCut = new TH1D("h_mmissCut","Missing mass with lead cuts;missing mass [GeV/c^{2}]",100,0,3);
  //plot the missing energy
  TH1D * h_emiss = new TH1D("h_emiss","Missing energy;missing energy",100,-1,2);
  TH1D * h_emissSmear = new TH1D("h_emissSmear","Missing energy;missing energy",100,-0.3,2);
  TH1D * h_emissCut = new TH1D("h_emissCut","Missing energy with lead cuts;missing energy",100,-1,2);

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
  

  //read in the data and make output file  
  TFile *f = new TFile("../../../data/1apr20/skim_C12_p.root");
  TFile *fout = new TFile("output/out_c12_p.root","RECREATE");

  //read in the Pmiss plots from eg2
  TFile *f2 = new TFile("../../../pmiss_eg2.root");
  TH1F *eg2_pmiss;
  TCanvas *c1 = (TCanvas*)f2->Get("Canvas_1");
  eg2_pmiss = (TH1F*)c1->GetPrimitive("h1");

  TCanvas *canvas = new TCanvas("canvas","output c12 p plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/out_c12_p.pdf";

  TTree * intree = (TTree*)f->Get("T");
  const int maxParticles=50;
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

  const int counts = 100000;
  double pmissO[counts];
  double pmissSmear[counts];
  double mmissO[counts];
  double mmissSmear[counts];
  int counter = 0;
  gRandom = new TRandom3();
  double p_over_q_O[counts];
  double p_over_q_sm[counts];
  double thetaPQ[counts];
  double thetaPQsm[counts];


  //where run table stored:
  //"../../calibration_data/run_table.dat"
  Fiducial fid_params(4461,2250,5996,"12C",true);
  double p_over_q_prev = 0;
  double theta_prev = 0;

  for (int event =0 ; event < intree->GetEntries() ; event++)
  //for (int event =0 ; event < 100000 ; event++)
    {
      if (event%10000==0)
      	cout << "Working on event " << event <<endl;

      intree->GetEvent(event);

       // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //cut on xB
      if (xB<1.1){continue;}

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at protons
	  if (type[part] != 2212)
	    continue;

	  //vertex cut-target dep.
	  if (vtxZ[part]<4 || vtxZ[part]>7){continue;}
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Ratio of p/q
	  double p_over_q = this_mom.Mag() / q.Mag();
	  //if ((p_over_q < 0.62)||(p_over_q > 0.96)) {continue;}

	  // Angle between p and q
	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
	  if (theta_pq>25.0){continue;}

	  //if (leadingID != -1){cout<<"already lead particle!"<<endl;}
	  
	  //first proton of interest in event
	  if (leadingID==-1){p_over_q_prev=p_over_q; theta_prev=theta_pq; leadingID=part;}

	  /////////////////////////////////////////////////////////////////////
	  //Here I selected the proton most likely to pass lead if more than 1:
	  //current p/q is larger and angle is smaller than previous
	  if ((p_over_q>p_over_q_prev) && (theta_pq<theta_prev)){leadingID=part;}

	  //one passes p/q cut
	  else if ((p_over_q<0.62 || p_over_q_prev<0.62) && (p_over_q>p_over_q_prev)){leadingID=part;}

	  //both pass p/q cut or neither pass p/q cut
	  else if (abs(p_over_q-0.62)/0.62*abs(theta_pq-25)/25.0 < abs(p_over_q_prev-0.62)/0.62*abs(theta_prev-25)/25.0)
	    {leadingID=part;}
	  	  
	  /*
	  //current p/q smaller and angle is smaller than previous
	  else if ((p_over_q<p_over_q_prev) && (theta_pq<theta_prev)){cout<<"p/q smaller, angle smaller"<<endl;}
	  //current p/q larger and angle is larger than previous
	  else if ((p_over_q>p_over_q_prev) && (theta_pq>theta_prev)){cout<<"p/q larger, angle larger"<<endl;}
	  //smaller, larger
	  else if ((p_over_q<p_over_q_prev) && (theta_pq>theta_prev)){cout<<"p/q smaller, angle larger"<<endl;}
	  */
	  /////////////////////////////////////////////////////////////////////	  
	}

      // If we have a proton passing the leading cuts, then fill the histogram
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,6.0,12.0);

	  
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());

	  double emiss =  Nu - (ep-Mp) - (Nu + M12C - ep - sqrt(pow(Nu + M12C - ep,2)- mom_miss.Mag2()));
	  double mmiss = sqrt(pow(Nu+2*Mp-ep,2.0) - (q - mom_lead).Mag2());
	  	  
	  //smear the proton momentum resolution
	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
	  TVector3 u1 = mom_lead.Unit();
	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());
	  
	  TVector3 mom_missSmear = mom_smear - q;
	  double epsmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = Nu - (epsmear-Mp) - (Nu + M12C - epsmear - sqrt(pow(Nu+ M12C-epsmear,2.0)-mom_missSmear.Mag2()));
	  double mmissSm = sqrt(pow(Nu+2*Mp-epsmear,2.0) - (q - mom_smear).Mag2());

	  TVector3 pECUVW(ec_u[leadingID],ec_v[leadingID],ec_w[leadingID]);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 

	  
	  hfid_pre->Fill(phi, theta);
	  //proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}
	  hec_pre->Fill(ec_x[leadingID],ec_y[leadingID]);

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}

	  hec_post->Fill(ec_x[leadingID],ec_y[leadingID]);	  
	  hfid_post->Fill(phi, theta);
	  h_zvtx_precut->Fill(vtxZ[leadingID]);

	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();
	  double theta_curr = mom_lead.Angle(q)*180.0/M_PI;
	  if (p_over_q>0.62&&p_over_q<1.1&&mom_miss.Mag()<1.0&&mmiss<1.1){
	      h_xb->Fill(xB);
	      h_nu->Fill(Nu);
	      h_etheta->Fill(mom_e.Theta()*180./M_PI);
	      h_emom->Fill(mom_e.Mag());
	      if (xB>=1.2 && mom_miss.Mag()>0.3&&p_over_q<0.96 ){h_q2->Fill(QSq);}
	      h_pmissCut->Fill(mom_miss.Mag());
	      h_mmissCut->Fill(mmiss);
	      h_pmissCut->Fill(mom_miss.Mag());
	      h_nmomCut->Fill(mom_lead.Mag());
	      h_nthetaCut->Fill(theta);
	      h_thetapqCut->Fill(mom_lead.Angle(q)*180./M_PI);
	      h_emissCut->Fill(emiss);
	      h_smearCut->Fill(mom_lead.Mag()-mom_smear.Mag());
	    }

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

	    pmissO[counter] = mom_miss.Mag();
	    pmissSmear[counter] = mom_missSmear.Mag();
	    mmissO[counter] = mmiss;
	    mmissSmear[counter] = mmissSm;
	    p_over_q_O[counter] = mom_lead.Mag() / q.Mag();
	    p_over_q_sm[counter] = mom_smear.Mag() / q.Mag();
	    thetaPQsm[counter] = mom_smear.Angle(q);
	    thetaPQ[counter] = mom_lead.Angle(q);

	    hist_pq_thetapq_smear->Fill(mom_smear.Mag()/q.Mag(),mom_smear.Angle(q)*180./M_PI);
	    h_smear->Fill(mom_lead.Mag()-mom_smear.Mag());
	    h_pmissSmear->Fill(mom_missSmear.Mag());
	    h_emissSmear->Fill(emissSmear);
	    h_mmissSmear->Fill(mmissSm);
	    h_nmomSmear->Fill(mom_smear.Mag());
	    h_emVpmSm->Fill(mom_missSmear.Mag(),emissSmear);
	    pvps->Fill(mom_miss.Mag(),mom_missSmear.Mag());

	    
	    counter++;
	  
	}//end leading
    }//end event

  double pmiss_cut[30];
  for (int ii=0; ii<30; ii++){
    pmiss_cut[ii] = (double)ii*2.0/100.0;
  }
  
  double mmiss_val[7] = {1.1,1.15,1.175,1.2,1.225,1.25,1.3};
  double falsePos[7][30];
  double falseNeg[7][30];
  double minimum = 9999;
  int mmissCut = 0;
  int pmissCut = 0;
  double falsePosSub[7][30];
  double falseNegSub[7][30];
  double falsePosOff[7][30];
  double falseNegOff[7][30];

  cout<<"pmiss "<<" mmiss "<<" false neg "<<" false pos "<<" false neg offset "<<" false pos offset "<<" false neg offSub "<<" false pos offSub "<<endl;
  
  //loop missing mass
  for (int ii=0;ii<7; ii++){
    for (int jj=0; jj<30; jj++){
      double fneg = 0;
      double fpos = 0;
      double den = 0;
      double fnegSub = 0;
      double fposSub = 0;
      double denSub = 0;

      for (int kk=0;kk<counter; kk++){
	///////////////////////////////////////////////////////////////////
	//Traditional Calculation//////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	//count smeared events that pass total (p/q, pmiss, mmiss, theta_pq)
	if (pmissSmear[kk]>pmiss_cut[jj] && pmissSmear[kk]<1.0 && mmissSmear[kk]<mmiss_val[ii] && p_over_q_sm[kk] >= 0.62 && p_over_q_sm[kk] <= 1.1){
	  den++;
	  //false positive: non-smeared fails, smeared passes
	  if(pmissO[kk]<0.3 || pmissO[kk]>1.0 || mmissO[kk]>1.1 || p_over_q_O[kk] < 0.62|| p_over_q_O[kk] > 1.1){
	    fpos++;
	  }
	  
	}
	//false negative: smeared fails, non-smeared passes
	if((pmissO[kk]>0.3 && pmissO[kk]<1.0 && mmissO[kk]<1.1 && p_over_q_O[kk] >= 0.62 && p_over_q_O[kk] <= 1.1)
	   &&  (pmissSmear[kk]<pmiss_cut[jj] || mmissSmear[kk]>mmiss_val[ii]|| pmissSmear[kk]>1.0 || p_over_q_sm[kk] < 0.62 || p_over_q_sm[kk] > 1.1)){
	  fneg++;
	}


	///////////////////////////////////////////////////////////////////
	//Offset Calculation///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	if (pmissO[kk]>pmiss_cut[jj] && pmissO[kk]<1.0 && mmissO[kk]<mmiss_val[ii] && p_over_q_O[kk] >= 0.62 && p_over_q_O[kk] <= 1.1){
	  denSub++;
	  //false positive: non-smeared fails, smeared passes
	  if(pmissO[kk]<0.3 || pmissO[kk]>1.0 || mmissO[kk]>1.1 || p_over_q_O[kk] < 0.62|| p_over_q_O[kk] > 1.1){
	    fposSub++;
	  }
	  
	}
	//false negative: smeared fails, non-smeared passes
	if((pmissO[kk]>0.3 && pmissO[kk]<1.0 && mmissO[kk]<1.1 && p_over_q_O[kk] >= 0.62 && p_over_q_O[kk] <= 1.1)
	   &&  (pmissO[kk]<pmiss_cut[jj] || mmissO[kk]>mmiss_val[ii]|| pmissO[kk]>1.0 || p_over_q_O[kk] < 0.62 || p_over_q_O[kk] > 1.1)){
	  fnegSub++;
	}
	
      }//end kk 

      
      //false +ive: events that should not pass do pass
      falsePos[ii][jj] = fpos*100.0/den;
      falsePosSub[ii][jj] = fpos*100.0/den - fposSub*100.0/den;//Sub;
      falsePosOff[ii][jj] = fposSub*100.0/den;//Sub;
      //false -ive: events that should pass do not pass
      falseNeg[ii][jj] = fneg*100.0/den;
      falseNegSub[ii][jj] = fneg*100.0/den - fnegSub*100.0/den;//Sub;
      falseNegOff[ii][jj] = fnegSub*100.0/den;//Sub;

      cout<<"  "<<pmiss_cut[jj]<<"  "<<mmiss_val[ii]<<"  "<<falseNeg[ii][jj]<<"  "<<falsePos[ii][jj]<<"  "<<falseNegOff[ii][jj]<<"  "<<falsePosOff[ii][jj]<<"  "<<falseNegSub[ii][jj]<<"  "<<falsePosSub[ii][jj]<<endl;

      if (abs(falsePos[ii][jj]-falseNeg[ii][jj])<minimum){
	minimum = abs(falsePos[ii][jj]-falseNeg[ii][jj]);
	mmissCut = ii;
	pmissCut = jj;
      }
      
    }//end jj
  }//end ii
  cout<<endl;
  cout<<"optimal pmiss cut: "<<pmiss_cut[pmissCut]<<" optimal mmiss cut: "<<mmiss_val[mmissCut]<<" false Pos: "<<falsePos[mmissCut][pmissCut]<<" false Neg: "<<falseNeg[mmissCut][pmissCut]<<" corrected false Pos: "<<falsePosSub[mmissCut][pmissCut]<<" corrected false Neg: "<<falseNegSub[mmissCut][pmissCut]<<" offset false Pos: "<<falsePosOff[mmissCut][pmissCut]<<" offset false Neg: "<<falseNegOff[mmissCut][pmissCut]<<endl;

  

  //normal false pos and false neg plot
  TGraph *grP[7];
  TGraph *grN[7];
  TMultiGraph *mgP = new TMultiGraph();
  TMultiGraph *mgN = new TMultiGraph();
  TLegend* legendP = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendN = new TLegend(0.75,0.6,0.85,0.8);
  legendP->SetHeader("m_{miss}");
  legendP->SetBorderSize(0);
  legendN->SetHeader("m_{miss}");
  legendN->SetBorderSize(0);
  mgP->SetTitle("False Positive;P_{miss} lower bound;False Positive [%]");
  mgN->SetTitle("False Negative;P_{miss} lower bound;False Negative [%]");
  for (int ii=0; ii<7; ii++){
    grP[ii] = new TGraph(30,pmiss_cut,falsePos[ii]);
    grP[ii]->SetMarkerStyle(20+ii);
    grP[ii]->SetMarkerColor(1+ii);
    mgP->Add(grP[ii],"p");
    legendP->AddEntry(grP[ii],Form("%.3f",mmiss_val[ii]),"p");
  
    grN[ii] = new TGraph(30,pmiss_cut,falseNeg[ii]);
    grN[ii]->SetMarkerStyle(20+ii);
    grN[ii]->SetMarkerColor(1+ii);
    mgN->Add(grN[ii],"p");
    legendN->AddEntry(grN[ii],Form("%.3f",mmiss_val[ii]),"p");
  }

  //offset false pos and false neg plot
  TGraph *grPoff[7];
  TGraph *grNoff[7];
  TMultiGraph *mgPoff = new TMultiGraph();
  TMultiGraph *mgNoff = new TMultiGraph();
  TLegend* legendPoff = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendNoff = new TLegend(0.75,0.6,0.85,0.8);
  legendPoff->SetHeader("m_{miss}");
  legendPoff->SetBorderSize(0);
  legendNoff->SetHeader("m_{miss}");
  legendNoff->SetBorderSize(0);
  mgPoff->SetTitle("False Positive, offset;P_{miss} lower bound;False Positive [%]");
  mgNoff->SetTitle("False Negative, offset;P_{miss} lower bound;False Negative [%]");
  for (int ii=0; ii<7; ii++){
    grPoff[ii] = new TGraph(30,pmiss_cut,falsePosOff[ii]);
    grPoff[ii]->SetMarkerStyle(20+ii);
    grPoff[ii]->SetMarkerColor(1+ii);
    mgPoff->Add(grPoff[ii],"p");
    legendPoff->AddEntry(grPoff[ii],Form("%.3f",mmiss_val[ii]),"p");
  
    grNoff[ii] = new TGraph(30,pmiss_cut,falseNegOff[ii]);
    grNoff[ii]->SetMarkerStyle(20+ii);
    grNoff[ii]->SetMarkerColor(1+ii);
    mgNoff->Add(grNoff[ii],"p");
    legendNoff->AddEntry(grNoff[ii],Form("%.3f",mmiss_val[ii]),"p");
  }


  //corrected false pos and false neg plot
  TGraph *grPsub[7];
  TGraph *grNsub[7];
  TMultiGraph *mgPsub = new TMultiGraph();
  TMultiGraph *mgNsub = new TMultiGraph();
  TLegend* legendPsub = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendNsub = new TLegend(0.75,0.6,0.85,0.8);
  legendPsub->SetHeader("m_{miss}");
  legendPsub->SetBorderSize(0);
  legendNsub->SetHeader("m_{miss}");
  legendNsub->SetBorderSize(0);
  mgPsub->SetTitle("False Positive, offset subtracted;P_{miss} lower bound;False Positive [%]");
  mgNsub->SetTitle("False Negative, offset subtracted;P_{miss} lower bound;False Negative [%]");
  for (int ii=0; ii<7; ii++){
    grPsub[ii] = new TGraph(30,pmiss_cut,falsePosSub[ii]);
    grPsub[ii]->SetMarkerStyle(20+ii);
    grPsub[ii]->SetMarkerColor(1+ii);
    mgPsub->Add(grPsub[ii],"p");
    legendPsub->AddEntry(grPsub[ii],Form("%.3f",mmiss_val[ii]),"p");
  
    grNsub[ii] = new TGraph(30,pmiss_cut,falseNegSub[ii]);
    grNsub[ii]->SetMarkerStyle(20+ii);
    grNsub[ii]->SetMarkerColor(1+ii);
    mgNsub->Add(grNsub[ii],"p");
    legendNsub->AddEntry(grNsub[ii],Form("%.3f",mmiss_val[ii]),"p");
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

  mgPoff->Draw("a");
  mgPoff->GetYaxis()->SetRangeUser(0,100);
  legendPoff->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgNoff->Draw("a");
  mgNoff->GetYaxis()->SetRangeUser(0,100);
  legendNoff->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
   mgPsub->Draw("a");
  mgPsub->GetYaxis()->SetRangeUser(-30,100);
  legendPsub->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgNsub->Draw("a");
  mgNsub->GetYaxis()->SetRangeUser(-30,100);
  legendNsub->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  pvps->Draw("colz");
  canvas->Print( (pdf_file_name + ")").c_str());



  fout->Write();
  fout->Close();
  cout<<"Good events: "<<totEvents<<endl;

}
