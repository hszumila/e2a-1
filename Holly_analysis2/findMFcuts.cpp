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
#include "Hconstants.h"
#include "Hfunctions.h"

//input is 3He, 4He, or 12C
int main(int argc, char ** argv){

  //////////////////////////////////
  // Set some target dep parameters
  /////////////////////////////////
  TString target = argv[1];
  cout<<"Target is: "<<target<<endl;

  //Change for fiducial cuts
  Fiducial fid_params(4461,2250,5996,"12C",true);

  //read in the data and make output file
  TFile *f = new TFile("../../../data/19jun20/skim_C12_p_Jun20.root");
  TFile *fout = new TFile("output_cuts/outMF_C12_p.root","RECREATE");
 
  std::string pdf_file_name = "output/outMF_C12_p.pdf";
   
  double tgt_min, tgt_max, tgt_M, tgt_R;
  int tgt_Z, tgt_A;
  if (target=="12C"){tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C; tgt_R = M11B;}
  else if (target=="3He"){tgt_min = vtx_min_He3; tgt_max = vtx_max_He3; tgt_Z = 2; tgt_A = 3; tgt_M = M3He; tgt_R = M2H;}
  else if (target=="4He"){tgt_min = vtx_min_He4; tgt_max = vtx_max_He4; tgt_Z = 2; tgt_A = 4; tgt_M = M4He; tgt_R = M3H;}

  cout<<"target min: "<<tgt_min<<" target max: "<<tgt_max<<endl;
  
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  //////////////////////////////////
  // Plots/////////////////////////
  /////////////////////////////////
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
  //plot y vs nu
  TH2F *h_yvnu_pre = new TH2F("h_yvnu_pre","before Em, Pm cuts;#nu;y",100,0,4.5,100,-1,1);
  TH2F *h_yvnu_post = new TH2F("h_yvnu_post","Em<140 MeV, Pm<220 MeV cuts;#nu;y",100,0,4.5,100,-1,1);
  //plot thetaPQ vs Q2
  TH2F *h_thetavq2_pre = new TH2F("h_thetavq2_pre","before Em, Pm cuts;Q^{2};#theta_{pq}",100,0,6,100,0,100);
  TH2F *h_thetavq2_post = new TH2F("h_thetavq2_post","Em<140 MeV, Pm<220 MeV cuts;Q^{2};#theta_{pq}",100,0,6,100,0,100);
  //plot x vs y
  TH2F *h_xvy_pre = new TH2F("h_xvy_pre","before Em, Pm cuts;x_{B};y",100,0,2,100,-1,1);
  TH2F *h_xvy_post = new TH2F("h_xvy_post","Em<140 MeV, Pm<220 MeV cuts;x_{B};y",100,0,2,100,-1,1);
  //plot x vs nu
  TH2F *h_xvnu_pre = new TH2F("h_xvnu_pre","before Em, Pm cuts;#nu;x_{B}",100,0,3,100,0,2);
  TH2F *h_xvnu_post = new TH2F("h_xvnu_post","Em<140 MeV, Pm<220 MeV cuts;#nu;x_{B}",100,0,3,100,0,2);
  //plot p vs psmear
  TH2F *pvps = new TH2F("h_pvps",";p_{miss} [GeV];p_{miss, from smeared}[GeV]",100,0,1,100,0,1);
  //plot nu vs Q2
  TH2F *h_nuVq_pre = new TH2F("h_nuVq_pre",";Q^{2};#nu",100,0,6,100,0,3);
  TH2F *h_nuVq_post = new TH2F("h_nuVq_post",";Q^{2};#nu",100,0,6,100,0,3);
  //plot Em vs Pm again for the selected QE cuts
  TH2F* h_emVpm_cuts =new TH2F("h_emVpm_cuts","found QE n cuts;pmiss;emiss",100,0,1,100,-0.8,0.8);
  TH2F* h_emVpm_sm_cuts =new TH2F("h_emVpm_sm_cuts","smeared,found QE n cuts;pmiss;emiss",100,0,1,100,-0.8,0.8);

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
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles], ec_in[nParticles], ec_out[nParticles];
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
  double emissO[counts];
  double emissSm[counts];
  int counter = 0;
  gRandom = new TRandom3();

  //where run table stored:
  //"../../calibration_data/run_table.dat"
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
      
      //mean field stuff:
      double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

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

	  // Angle between p and q
	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
	  if (theta_pq>25.0){continue;}
	  
	  //first proton of interest in event
	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

	  /////////////////////////////////////////////////////////////////////
	  //Here I selected the proton with angle closest to q-vector:	  
	  //current p/q is larger and angle is smaller than previous
	  if (theta_pq<theta_prev){leadingID=part;}
	  /////////////////////////////////////////////////////////////////////	  
	 
	}

      // If we have a proton passing the leading cuts, then fill the histogram
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  if (mom_lead.Mag()>1.95){continue;}
	  
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

	  //if (mom_smear.Mag()>2.34){continue;}
	  
	  TVector3 mom_missSmear = mom_smear - q;
	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);

	  //TVector3 pECUVW(ec_u[leadingID],ec_v[leadingID],ec_w[leadingID]);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 

	  //cut at Pmiss>1:
	  //if (mom_miss.Mag()>1.0 || mom_missSmear.Mag()>1.0){continue;}
	  
	  hfid_pre->Fill(phi, theta);
	  //proton fiducial cuts:
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

	 
	  ////////////////////////////////////////////////
	  //plots for determining electron cutsfor QE:
	  h_yvnu_pre->Fill(Nu,y_sc);
	  h_thetavq2_pre->Fill(QSq,mom_lead.Angle(q)*180./M_PI);
	  h_xvy_pre->Fill(xB,y_sc);
	  h_xvnu_pre->Fill(Nu,xB);
	  h_nuVq_pre->Fill(QSq,Nu);
	  hist_pq_thetapq->Fill(mom_lead.Mag()/q.Mag(),mom_lead.Angle(q)*180./M_PI);
	  hist_pq_thetapq_smear->Fill(mom_smear.Mag()/q.Mag(),mom_smear.Angle(q)*180./M_PI);
	  hist_pmiss_xb->Fill(mom_miss.Mag(),xB);

	  
	  if ( mom_miss.Mag()<pmiss_MFnomCut && emiss<emiss_nomCut){
	    h_yvnu_post->Fill(Nu,y_sc);
	    h_thetavq2_post->Fill(QSq,mom_lead.Angle(q)*180./M_PI);
	    h_xvy_post->Fill(xB,y_sc);
	    h_xvnu_post->Fill(Nu,xB);
	    h_nuVq_post->Fill(QSq,Nu);
	  }
	  ///////////////////////////////////////////////

	  //plot the (e,e'p) MF selected events and count them
	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut && mom_miss.Mag()<pmiss_MFnomCut && emiss<emiss_nomCut){	       
	    h_xb->Fill(xB);
	    h_nu->Fill(Nu);
	    h_etheta->Fill(mom_e.Theta()*180./M_PI);
	    h_emom->Fill(mom_e.Mag());
	    h_q2->Fill(QSq);
	    h_mmissCut->Fill(mmiss);
	    h_nmomCut->Fill(mom_lead.Mag());
	    h_nthetaCut->Fill(theta);
	    h_thetapqCut->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_emissCut->Fill(emiss);
	    h_smearCut->Fill(mom_lead.Mag()-mom_smear.Mag());	   
	  }
	  	  
	  //plot with all cuts except Em, Pm
	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut && counter<counts){
	    pvps->Fill(mom_miss.Mag(),mom_missSmear.Mag());
	    h_pmissCut->Fill(mom_miss.Mag());
	    h_mmissSmear->Fill(mmissSmear);
	    h_nmomSmear->Fill(mom_smear.Mag());
	    h_pmiss->Fill(mom_miss.Mag());
	    h_pmissSmear->Fill(mom_missSmear.Mag());
	    h_emissSmear->Fill(emissSmear);
	    h_emiss->Fill(emiss);
	    if (emiss<emiss_nomCut &&mom_miss.Mag()<pmiss_MFnomCut){
	      h_emVpmSm->Fill(mom_missSmear.Mag(),emissSmear);
	    }
	    h_emVpm->Fill(mom_miss.Mag(),emiss);
	    h_pmMpm->Fill(mom_lead.Mag()-mom_miss.Mag());
	    h_mmiss->Fill(mmiss);
	    h_nmom->Fill(mom_lead.Mag());
	    h_ntheta->Fill(theta);
	    h_thetapq->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_timeres->Fill(sc_time[leadingID]-ec_time[leadingID]-(sc_path[leadingID]-ec_path[leadingID])/cc);
	    h_smear->Fill(mom_lead.Mag()-mom_smear.Mag());
	      
	    pmissO[counter] = mom_miss.Mag();
	    pmissSm[counter] = mom_missSmear.Mag();
	    emissO[counter] = emiss;
	    emissSm[counter] = emissSmear;
	    counter++;
	  }
	  if(mom_missSmear.Mag()<0.28 && emissSmear<0.19){
	    if(mom_miss.Mag()>0.25 || emiss>0.08){
	      h_emVpm_cuts->Fill(mom_miss.Mag(),emiss);
	    }
	    h_emVpm_sm_cuts->Fill(mom_missSmear.Mag(),emissSmear);
	  }
	
	}//end leading
    }//end event
   
  //Calculate the false positives and false negatives.
  const int n_pmiss_cut = 30;
  const int n_emiss_cut =11;
  //make the Pmiss cuts:
  double pmiss_cut[n_pmiss_cut];
  for (int ii=0; ii<n_pmiss_cut; ii++){
    pmiss_cut[ii] = (double)ii*1.5/100.0+0.1;
  }

  //make Emiss cuts: 0.14 emiss, 0.22 pmiss
  //-used:double emiss_val[n_emiss_cut] = {0.18,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4};//{0.07,0.1,0.13,0.16,0.17,0.18,0.19,0.20,0.22,0.24,0.26};
  double emiss_val[n_emiss_cut] = {0.1,0.15,0.20,0.23,0.25,0.265,0.28,0.295,0.31,0.34,0.37};
  //double emiss_val[n_emiss_cut] = {0.13,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.3,0.32};//{0.07,0.1,0.13,0.16,0.17,0.18,0.19,0.20,0.22,0.24,0.26};
  double falsePos[n_emiss_cut][n_pmiss_cut];
  double falseNeg[n_emiss_cut][n_pmiss_cut];  
  double minimum = 9999;
  int emissCut = 0;
  int pmissCut = 0;
  double falsePosO[n_emiss_cut][n_pmiss_cut];
  double falseNegO[n_emiss_cut][n_pmiss_cut];
  double falsePosErr[n_emiss_cut][n_pmiss_cut];
  double falseNegErr[n_emiss_cut][n_pmiss_cut];

 cout<<"pmiss: "<<" emiss: "<<" false neg: "<<" false pos: "<<" false neg error: "<<" false pos error: "<<" fneg-fpos: "<<" err(fneg-fpos): "<<endl;
 
  //loop missing energy
  for (int ii=0;ii<n_emiss_cut; ii++){
    for (int jj=0; jj<n_pmiss_cut; jj++){
      double fneg = 0;
      double fpos = 0;
      double den = 0;
      double fnegO = 0;
      double fposO = 0;
      double fnegSub = 0;
      double fposSub = 0;
      double truNeg = 0;
      double truPos = 0;

      for (int kk=0;kk<counter; kk++){
	///////////////////////////////////////////////////////////////////
	//Traditional Calculation//////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	//count smeared events that pass total (pmiss, emiss, theta_pq)
	if (pmissSm[kk]<pmiss_cut[jj] && emissSm[kk]<emiss_val[ii]){
	  den++;
	  //false positive: non-smeared fails, smeared passes
	  if( pmissO[kk]>pmiss_MFnomCut || emissO[kk]>emiss_nomCut ){
	    fposO++;
	  }
	  
	}
	//false negative: smeared fails, non-smeared passes
	if((pmissO[kk]<pmiss_MFnomCut && emissO[kk]<emiss_nomCut)
	  &&  (pmissSm[kk]>pmiss_cut[jj] || emissSm[kk]>emiss_val[ii])){
	  fnegO++;
	}

	///////////////////////////////////////////////////////////////////	
	//false pos
	if ((pmissSm[kk]<pmiss_cut[jj] && emissSm[kk]<emiss_val[ii])
	    &&(pmissO[kk]>pmiss_MFnomCut || emissO[kk]>emiss_nomCut)){
	  fpos++;
	}
	//true negative
	if (pmissO[kk]>pmiss_MFnomCut || emissO[kk]>emiss_nomCut){
	  truNeg++;
	}
	//false neg
	if((pmissO[kk]<pmiss_MFnomCut && emissO[kk]<emiss_nomCut)
	  &&  (pmissSm[kk]>pmiss_cut[jj] || emissSm[kk]>emiss_val[ii])){
	  fneg++;
	}
	//true positive
	if (pmissO[kk]<pmiss_MFnomCut && emissO[kk]<emiss_nomCut){
	  truPos++;
	}	
      }//end kk 
      
      falsePos[ii][jj] = fpos*100.0/(fpos+truNeg);
      falsePosO[ii][jj] = fposO*100.0/den;
      falsePosErr[ii][jj] = sqrt(fposO+pow(falsePosO[ii][jj]/100.0,2.0)*den)/den;

      falseNeg[ii][jj] = fneg*100/(truPos+fneg);
      falseNegO[ii][jj] = fnegO*100/den;
      falseNegErr[ii][jj] = sqrt(fnegO+pow(falseNegO[ii][jj]/100.0,2.0)*den)/den;
      double diff = falseNegO[ii][jj] - falsePosO[ii][jj];
      double errDiff = sqrt(pow(falsePosErr[ii][jj],2.0)+pow(falseNegErr[ii][jj],2.0))*100.0;
      
        cout<<pmiss_cut[jj]<<","<<emiss_val[ii]<<","<<falseNegO[ii][jj]<<","<<falsePosO[ii][jj]<<","<<falseNegErr[ii][jj]<<","<<falsePosErr[ii][jj]<<","<<diff<<","<<errDiff<<endl;
      
      if (abs(falsePosO[ii][jj]-falseNegO[ii][jj])<minimum){
	minimum = abs(falsePosO[ii][jj]-falseNegO[ii][jj]);
	emissCut = ii;
	pmissCut = jj;
      }
      
    }//end jj
  }//end ii
  cout<<"optimal pmiss cut: "<<pmiss_cut[pmissCut]<<" optimal emiss cut: "<<emiss_val[emissCut]<<" false Pos: "<<falsePosO[emissCut][pmissCut]<<" false Neg: "<<falseNegO[emissCut][pmissCut]<<endl;
  

  //normal false pos and false neg plot
  TGraph *grP[n_emiss_cut];
  TGraph *grN[n_emiss_cut];
  TMultiGraph *mgP = new TMultiGraph();
  TMultiGraph *mgN = new TMultiGraph();
  TLegend* legendP = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendN = new TLegend(0.75,0.6,0.85,0.8);
  legendP->SetHeader("E_{miss}");
  legendP->SetBorderSize(0);
  legendN->SetHeader("E_{miss}");
  legendN->SetBorderSize(0);
  mgP->SetTitle("False Positive;P_{miss} upper cut;False Positive [%]");
  mgN->SetTitle("False Negative;P_{miss} upper cut;False Negative [%]");
  for (int ii=0; ii<n_emiss_cut; ii++){
    grP[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falsePos[ii]);
    grP[ii]->SetMarkerStyle(20+ii);
    grP[ii]->SetMarkerColor(1+ii);
    mgP->Add(grP[ii],"p");
    legendP->AddEntry(grP[ii],Form("%.3f",emiss_val[ii]),"p");
  
    grN[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falseNeg[ii]);
    grN[ii]->SetMarkerStyle(20+ii);
    grN[ii]->SetMarkerColor(1+ii);
    mgN->Add(grN[ii],"p");
    legendN->AddEntry(grN[ii],Form("%.3f",emiss_val[ii]),"p");
  }

  //false pos and false neg plot with smeared passing denominator
  TGraph *grPO[n_emiss_cut];
  TGraph *grNO[n_emiss_cut];
  TMultiGraph *mgPO = new TMultiGraph();
  TMultiGraph *mgNO = new TMultiGraph();
  TLegend* legendPO = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendNO = new TLegend(0.75,0.6,0.85,0.8);
  legendPO->SetHeader("E_{miss}");
  legendPO->SetBorderSize(0);
  legendNO->SetHeader("E_{miss}");
  legendNO->SetBorderSize(0);
  mgPO->SetTitle("False Positive Rate;P_{miss} upper cut;False Positive Rate [%]");
  mgNO->SetTitle("False Negative per smeared passing;P_{miss} upper cut;False Negative per smeared passing [%]");
  for (int ii=0; ii<n_emiss_cut; ii++){
    grPO[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falsePosO[ii]);
    grPO[ii]->SetMarkerStyle(20+ii);
    grPO[ii]->SetMarkerColor(1+ii);
    mgPO->Add(grPO[ii],"p");
    legendPO->AddEntry(grPO[ii],Form("%.3f",emiss_val[ii]),"p");
  
    grNO[ii] = new TGraph(n_pmiss_cut,pmiss_cut,falseNegO[ii]);
    grNO[ii]->SetMarkerStyle(20+ii);
    grNO[ii]->SetMarkerColor(1+ii);
    mgNO->Add(grNO[ii],"p");
    legendNO->AddEntry(grNO[ii],Form("%.3f",emiss_val[ii]),"p");
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
  
  h_yvnu_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_yvnu_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetavq2_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetavq2_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xvy_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xvy_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xvnu_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xvnu_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  pvps->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nuVq_pre->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nuVq_post->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm_cuts->Draw("colz");
  cout<<"Non-smear bad: "<<h_emVpm_cuts->Integral()<<endl;
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm_sm_cuts->Draw("colz");
  cout<<"Smeared good: "<<h_emVpm_sm_cuts->Integral()<<endl;
  canvas->Print( (pdf_file_name + ")").c_str());


  fout->Write();
  fout->Close();
  cout<<"Good events: "<<totEvents<<endl;

}
