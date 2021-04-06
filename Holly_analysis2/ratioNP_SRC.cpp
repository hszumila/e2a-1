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
#include "Acceptance.h"
#include "Hconstants.h"
#include "Hfunctions.h"


int main(){

  int event_p = 0;
  int event_n = 0;
  int stats_p = 0;
  int stats_n = 0;
  double tgt_min, tgt_max, tgt_M, tgt_R;
  int tgt_Z, tgt_A;
  const int maxParticles=50;
  
  //input 12C
  TString target = "12C";
  if (target=="12C"){tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C;}
  else if (target=="3He"){tgt_min = vtx_min_He3; tgt_max = vtx_max_He3; tgt_Z = 2; tgt_A = 3; tgt_M = M3He; tgt_R = M2H;}
  else if (target=="4He"){tgt_min = vtx_min_He4; tgt_max = vtx_max_He4; tgt_Z = 2; tgt_A = 4; tgt_M = M4He; tgt_R = M3H;}
  Fiducial fid_params(4461,2250,5996,"12C",true);
  Acceptance myMap("C12", 4461, 2250, "p");

  //neutron input file
  TFile *f_3n = new TFile("../../../data/19jun20/skim_C12_n_Mar5.root");  
  //TFile *f_3n = new TFile("../../../data/p3/skim_pass3_n.root");  
  //proton input file
  TFile *f_3p = new TFile("../../../data/19jun20/skim_C12_p_Jun20.root");
  //TFile *f_3p = new TFile("../../../data/p3/skim_pass3_p.root");

  
  //output file and plots
  TString outputpdf = "output/output_npratio_SRC_12C";

    
  //plot the proton momentum with XS weighting
  TH1F* h_p_XS =new TH1F("h_p_XS","proton momentum; proton momentum [GeV/c];1/#sigma_{p}",10,0,2.34);
  TH1F* h_psm_XS =new TH1F("h_psm_XS","smeared proton momentum; smeared proton momentum [GeV/c];1/#sigma_{p}",10,0,2.34);
    
  //plot the neutron momentum with XS weighting
  TH1F* h_n_XS =new TH1F("h_n_XS","neutron momentum; neutron momentum [GeV/c];1/#sigma_{n}",10,0,2.34);
  
  //plot the proton momentum without XS weighting
  TH1F* h_p =new TH1F("h_p","proton momentum; proton momentum [GeV/c]",10,0,2.34);
  TH1F* h_psm =new TH1F("h_psm","smeared proton momentum; smeared proton momentum [GeV/c]",10,0,2.34);

  //plot the neutron momentum without XS weighting
  TH1F* h_n =new TH1F("h_n","neutron momentum; neutron momentum [GeV/c]",10,0,2.34);

  //plot the n/p ratio
  TH1F *h_ratio = new TH1F("h_ratio","n/p; momentum [GeV/c]; n/p",10,0,2.34);
  
  //plot proton acceptance
  TH2F *h2_p = new TH2F("h2_p","proton;#phi [deg];#theta [deg]",100,-100,360,100,0,70);

  //plot neutron acceptance
  TH2F *h2_n = new TH2F("h2_n","neutron;#phi [deg];#theta [deg]",100,-100,360,100,0,70);
  
  TTree * intree = (TTree*)f_3p->Get("T");
  int nParticles, type[maxParticles],nRun;
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles],ec_in[maxParticles],ec_out[maxParticles];
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
  intree->SetBranchAddress("nRun",&nRun);

  double theta_prev = 0;
  double p_over_q_prev = 0;

  //get (e,e'p)
  for (int event =0 ; event < intree->GetEntries() ; event++)
    //for (int event =0 ; event < 100000 ; event++)
    {
      if (event%10000==0){cout << "Working on event " << event <<endl;}
      
      intree->GetEvent(event);
      //if (nRun!=18020){continue;}

      // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //cut on xB
      if (xB<xb_cut){continue;}

      if (QSq<Q2_min){continue;}
      
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
	  
	}
      
      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);
	  
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
	  double p_over_q_sm = mom_smear.Mag()/q.Mag();
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_smear)) {continue;}
	  
	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_smear,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_smear.Mag()>1.95) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_smear.Mag() / q.Mag();

	  //check if in good acceptance
	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}
	  
	  //select the 12C proton sample
	  if (p_over_q_sm>pq_min && p_over_q_sm<pq_max_n 
			  && mmissSmear<mmiss_nomCut_n && mom_missSmear.Mag()<1 && mom_missSmear.Mag()>pmissSRC_nomCut_n
	  //if (p_over_q>pq_min && p_over_q<pq_max_n 
		//	  && mmiss<mmiss_nomCut && mom_miss.Mag()<1 && mom_miss.Mag()>pmiss_nomCut
			  //&& mom_smear.Mag()>1.4 && mom_smear.Mag()<1.8
			  ){
	    
	    double weight = getSigmaWeight(mom_e.Theta(), QSq, 1);
	   
	    double p_eff = myMap.get_acc(mom_lead);

	    event_p += 1./p_eff * 1./weight;
	    stats_p += 1;
	    h_p_XS->Fill(mom_lead.Mag(),1./p_eff * 1./weight);
	    h_psm_XS->Fill(mom_smear.Mag(),1./p_eff * 1./weight);
	    h_p->Fill(mom_lead.Mag());
	    h_psm->Fill(mom_smear.Mag());
	    h2_p->Fill(phi,theta);

	    
	  }

	}//end leading
    }//end event
  
  intree->Delete();
  
  //get (e,e'n)
  //open the neutron file
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
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  intree->SetBranchAddress("nRun",&nRun);
  
  for (int event =0 ; event < intree->GetEntries() ; event++)
    //for (int event =0 ; event < 100000 ; event++)
    {
      if (event%10000==0){cout << "Working on event " << event <<endl;}

      intree->GetEvent(event);
      //if (nRun!=18020){continue;}

       // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //cut on xB
      if (xB<xb_cut){continue;}

      if (QSq<Q2_min){continue;}

      // Loop over particles looking for a neutron  most likely to pass lead cuts
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
	  TVector3 nmomCorr = correctPnVec(this_mom, ec_in[part], ec_out[part]);
	  double p_over_q = nmomCorr.Mag() / q.Mag();
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
	  	  
	}

      //Now we have a leading neutron!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  mom_lead = correctPnVec(mom_lead,ec_in[leadingID],ec_out[leadingID]);
	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

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
	  if(mom_lead.Mag()>1.95) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();

	  //check if in good acceptance
	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}

	  //select the 3He neutron sample
	  if (p_over_q>pq_min && p_over_q<pq_max_n && mmiss<mmiss_nomCut_n 
			  && mom_miss.Mag()<1 && mom_miss.Mag()>pmissSRC_nomCut_n
			  //&& mom_lead.Mag()>1.4 && mom_lead.Mag()<1.8
			  ){

	    double weight = getSigmaWeight(mom_e.Theta(), QSq, 2);	    
	    double n_eff = getNEff(mom_lead.Mag());

	    event_n += 1./(n_eff) * 1./weight;
	    stats_n += 1;
	    
	    h_n_XS->Fill(mom_lead.Mag(),1./n_eff * 1./weight);
	   	h_ratio->Fill(mom_lead.Mag(),1./n_eff * 1./weight);
	    h_n->Fill(mom_lead.Mag());
	    h2_n->Fill(phi,theta);


	  }	  

	}//end leading
    }//end event


  double ratio = (double) event_n/event_p;
  double stats = sqrt(pow(1./stats_n,2.0)+pow(1./stats_p,2.0));
  cout<<"proton: "<<event_p<<" neutron: "<<event_n<<" ratio: "<<ratio<<" stats p: "<<stats_p<<" stats n: "<<stats_n<<endl;

  TCanvas * c1 = new TCanvas("c1","",1200,900);
   c1->Divide(3,2);
   c1->cd(1); h_p->Draw();
   c1->cd(2); h_psm->Draw();
   c1->cd(3); h_n->Draw();
   c1->cd(4); h_p_XS->Draw("hist");
   c1->cd(5); h_psm_XS->Draw("hist");
   c1->cd(6); h_n_XS->Draw("hist");
   c1->Print(outputpdf+".pdf(");

   //plot the ratio n/p vs momentum 
   TCanvas * c2 = new TCanvas("c2","",1200,900);
   h_ratio->Divide(h_psm_XS);
   /*
     for( int i = 1 ; i <= h_psm_XS -> GetXaxis()->GetNbins(); i++){
   	  double A  = h_psm_XS -> GetBinContent(i);
   	  double B  = h_n_XS -> GetBinContent(i);
   	  double eA = h_psm_XS -> GetBinError  (i);
   	  double eB = h_n_XS -> GetBinError  (i);
   	  double error_both = sqrt( pow( eA/B ,2) + pow( A*eB/ pow( B  ,2)  ,2) );
   	  h_ratio -> SetBinError(i,error_both);
     }
     */
     c2->cd(); h_ratio->Draw("hist");
     c2->Print(outputpdf+".pdf(");
     
     
     TCanvas * c3 = new TCanvas("c3","",1200,900);
     c3->Divide(2,1);
     c3->cd(1); h2_p->Draw("colz");
     c3->cd(2); h2_n->Draw("colz");
     c3->Print(outputpdf+".pdf)");

}
