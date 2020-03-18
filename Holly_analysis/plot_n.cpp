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

///////////////////////////////
//TODO:
//3. correct smear function
//////////////////////////////

int main(){

  gROOT->SetBatch(true);

  //plot theta vs phi before
  TH2F *hfid_pre = new TH2F("hfid_pre","no p fid cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot theta vs phi after
  TH2F *hfid_post = new TH2F("hfid_post","with p fid cut;#phi [deg];#theta [deg]",200,-100,360,200,0,70);
  //plot EC x vs y before
  TH2F *hec_pre = new TH2F("hec_pre","no fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
  //plot EC x vs y after
  TH2F *hec_post = new TH2F("hec_post","10cm fid cut; EC_{x} [cm];EC_{y} [cm]",200,-500,500,200,-500,500);
  //plot pmiss and xb
  TH2D * hist_pmiss_xb = new TH2D("pmiss_xb","Events passing leading cuts;pmiss;Xb;Counts",28,0.3,1.0,40,1,3);
  //plot theta_pq and p/q
  TH2D * hist_pq_thetapq = new TH2D("hist_pq_thetapq","Events passing leading cuts;p/q;#theta_{pq} [deg];Counts",100,0,1.5,100,0,70);
  //plot z vertex
  TH1D * h_zvtx_precut = new TH1D("h_zvtx_precut","Corrected z vertex-no cut;z vertex [cm]",100,-10,10);
  //plot the missing mass
  TH1D * h_mmiss = new TH1D("h_mmiss","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  TH1D * h_mmissSmear = new TH1D("h_mmissSmear","Missing mass;missing mass [GeV/c^{2}]",100,0,3);
  //plot the missing energy
  TH1D * h_emiss = new TH1D("h_emiss","Missing energy;missing energy",100,-1,2);
  TH1D * h_emissSmear = new TH1D("h_emissSmear","Missing energy;missing energy",100,-0.3,2);
  //plot Q^2
  TH1F* h_q2 =new TH1F("h_q2",";Q^{2};",100,0,6);
  //plot xB
  TH1F* h_xb =new TH1F("h_xb",";x_{B};",100,0,2);
  //plot nu
  TH1F* h_nu =new TH1F("h_nu",";#nu [GeV];",100,0,3);
  //plot e- theta
  TH1F* h_etheta =new TH1F("h_etheta",";#theta_{e-} [deg];",100,0,45);
  //plot e- momentum
  TH1F* h_emom =new TH1F("h_emom",";P_{e-} [GeV/c];",100,0,6);
  //plot pmiss
  TH1F* h_pmiss =new TH1F("h_pmiss",";P_{miss} [GeV/c];",100,0,1.5);
  TH1F* h_pmissSmear =new TH1F("h_pmissSmear",";P_{miss} [GeV/c];",100,0,1.5);
  //plot nuclean momentum
  TH1F* h_nmom =new TH1F("h_nmom",";P_{lead N} [GeV/c];",100,0,3);
  TH1F* h_nmomSmear =new TH1F("h_nmomSmear",";P_{lead N} [GeV/c];",100,0,3);
  //plot theta_n
  TH1F* h_ntheta =new TH1F("h_ntheta",";#theta_{lead N} [deg];",100,0,90);
  //plot theta_pq
  TH1F* h_thetapq =new TH1F("h_thetapq",";#theta_{lead N/q} [deg];",100,0,25);
  //plot proton smearing
  TH1F* h_smear =new TH1F("h_smear",";p_{p} - p_{p,smeared} [GeV/c];",100,-1,1);
  //plot nominal quantities
  TH2F* h_emVpm =new TH2F("h_emVpm",";pmiss;emiss",100,0,0.5,100,0.0,0.5);
  //plot smeared quantities
  TH2F* h_emVpmSm =new TH2F("h_emVpmSm","smeared;pmiss;emiss",100,0,0.5,100,0.0,0.5);


  //read in the data and make output file  
  TFile *f = new TFile("../../../data/nina/skim_C12_n.root");
  TFile *fout = new TFile("output/out_c12_n.root","RECREATE");

  TCanvas *canvas = new TCanvas("canvas","output c12 n plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/out_c12_n.pdf";

  TTree * intree = (TTree*)f->Get("T");
  const int maxParticles=50;
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles];
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

  const int counts = 100000;
  double pmissO[counts];
  double pmissSmear[counts];
  double mmissO[counts];
  double mmissSmear[counts];
  int counter = 0;
  gRandom = new TRandom3();
  double p_over_q_O[counts];
  double p_over_q_sm[counts];


  //where run table stored:
  //"../../calibration_data/run_table.dat"
  cout << "Before fiducial class" << endl;
  Fiducial fid_params(4461,2250,5996,"12C",true);
  cout << "After fiducial class" << endl;
  double p_over_q_prev = 0;

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

      // Loop over particles looking for a proton that will pass leading cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at protons
	  if (type[part] != 2212)
	    continue;
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Cut on the angle between p and q
	  if (this_mom.Angle(q) > 25.*M_PI/180.) {continue;}

	  // Cut on the ratio of p/q
	  double p_over_q = this_mom.Mag() / q.Mag();
	  //if ((p_over_q < 0.62)||(p_over_q > 0.96)) {continue;}
	  //if ((p_over_q < 0.62)||(p_over_q > 1.1)) {continue;}
	  // This proton passes the cuts, it is the leading proton for this event
	  if (leadingID==-1){p_over_q_prev = p_over_q;}

	  //if (leadingID>0) {cout<<"Already have lead p!!!"<<endl;
	    //cout<<"     p/q current: "<<p_over_q<<endl;
	    //cout<<"     p/q previous: "<<p_over_q_prev<<endl;

	  if (p_over_q>p_over_q_prev){leadingID=part;}
	  else{continue;}
	  
	  }

	  
	  leadingID = part;
	}

      // If we have a proton passing the leading cuts, then fill the histogram
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());

	  //emiss and mmiss distr	  
	  //assumes ground state:
	  //double emiss = Nu + 12.0107*amu - ep - sqrt(pow(10.811*amu,2.0)+(q-mom_lead).Mag2());

	  //allows for excited state:
	  //double mmiss = sqrt(pow(emiss,2.0)-mom_miss.Mag2());

	  double emiss =  Nu - (ep-Mp) - (Nu + M12C - ep - sqrt(pow(Nu + M12C - ep,2)- mom_miss.Mag2()));
	  double mmiss = sqrt(pow(Nu+2*Mp-ep,2.0) - (q - mom_lead).Mag2());
	  
	  
	  //smear the proton momentum resolution
	  TVector3 mom_smear(mom_lead.X()+gRandom->Gaus(0.0, 0.11)*mom_lead.X(),
			     mom_lead.Y()+gRandom->Gaus(0.0, 0.11)*mom_lead.Y(),
			     mom_lead.Z()+gRandom->Gaus(0.0, 0.11)*mom_lead.Z());
	  TVector3 mom_missSmear = mom_smear - q;
	  double epsmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
	  double emissSmear = Nu - (epsmear-Mp) - (Nu + M12C - epsmear - sqrt(pow(Nu+ M12C-epsmear,2.0)-mom_missSmear.Mag2()));
	  double mmissSm = sqrt(pow(Nu+2*Mp-epsmear,2.0) - (q - mom_smear).Mag2());

	  TVector3 pECUVW(ec_u[leadingID],ec_v[leadingID],ec_w[leadingID]);

	  //cout<<"emiss: "<<emiss<<" emiss smeared: "<<emissSmear<<" mmiss: "<<mmiss<<" mmiss smeared: "<<mmissSm<<endl;
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 

	  /////////////////////////////////////////
	  //need to implement proton fiducial cuts
	  /////////////////////////////////////////
	  
	  hfid_pre->Fill(phi, theta);
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}
	  hec_pre->Fill(ec_x[leadingID],ec_y[leadingID]);	  
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  hec_post->Fill(ec_x[leadingID],ec_y[leadingID]);	  
	  hfid_post->Fill(phi, theta);
	  h_zvtx_precut->Fill(vtxZ[leadingID]);

	  if (vtxZ[leadingID]>4 && vtxZ[leadingID]<7){
	    h_q2->Fill(QSq);
	    hist_pmiss_xb->Fill(mom_miss.Mag(),xB);
	    h_mmiss->Fill(mmiss);
	    h_xb->Fill(xB);
	    h_nu->Fill(Nu);
	    h_etheta->Fill(mom_e.Theta()*180./M_PI);
	    h_emom->Fill(mom_e.Mag());
	    h_pmiss->Fill(mom_miss.Mag());
	    h_nmom->Fill(mom_lead.Mag());
	    h_ntheta->Fill(theta);
	    h_thetapq->Fill(mom_lead.Angle(q)*180./M_PI);
	    h_emiss->Fill(emiss);
	    hist_pq_thetapq->Fill(mom_lead.Mag()/q.Mag(),mom_lead.Angle(q)*180./M_PI);
	    h_emVpm->Fill(mom_miss.Mag(),emiss);

	    pmissO[counter] = mom_miss.Mag();
	    pmissSmear[counter] = mom_missSmear.Mag();
	    mmissO[counter] = mmiss;
	    mmissSmear[counter] = mmissSm;
	    p_over_q_O[counter] = mom_lead.Mag() / q.Mag();
	    p_over_q_sm[counter] = mom_smear.Mag() / q.Mag();

	    h_smear->Fill(mom_lead.Mag()-mom_smear.Mag());
	    h_pmissSmear->Fill(mom_missSmear.Mag());
	    h_emissSmear->Fill(emissSmear);
	    h_mmissSmear->Fill(mmissSm);
	    h_nmomSmear->Fill(mom_smear.Mag());
	    h_emVpmSm->Fill(mom_missSmear.Mag(),emissSmear);

	    
	    counter++;
	  }
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
  
  //loop missing mass
  for (int ii=0;ii<7; ii++){
    for (int jj=0; jj<30; jj++){
      double fneg = 0;
      double fpos = 0;
      double den = 0;
      
      for (int kk=0;kk<counter; kk++){
	//count smeared events that pass total
	if (pmissSmear[kk]>pmiss_cut[jj] && mmissSmear[kk]<mmiss_val[ii] && p_over_q_sm[kk] >= 0.62 && p_over_q_sm[kk] <= 0.96){
	  den++;
	  //false positive: non-smeared fails, smeared passes
	  if(pmissO[kk]<0.3 || mmissO[kk]>1.1 || p_over_q_O[kk] < 0.62|| p_over_q_O[kk] > 0.96){
	    fpos++;
	  }
	  
	}
	//false negative: smeared fails, non-smeared passes
	  if((pmissO[kk]>0.3 && mmissO[kk]<1.1 && p_over_q_O[kk] >= 0.62 && p_over_q_O[kk] <= 0.96) &&  (pmissSmear[kk]<pmiss_cut[jj] || mmissSmear[kk]>mmiss_val[ii])){
	  fneg++;
	}		
      }//end kk 
          
      //false +ive: events that should not pass do pass
      falsePos[ii][jj] = fpos*100.0/den;
      //false -ive: events that should pass do not pass
      falseNeg[ii][jj] = fneg*100.0/den;
      //cout<<"pmiss low cut: "<<pmiss_cut[jj]<<" total smeared passing: "<<den<<" false positives: "<<fpos<<" false negatives: "<<fneg<<" all uncut events: "<<counter<<endl;
      if (abs(falsePos[ii][jj]-falseNeg[ii][jj])<minimum){
	minimum = abs(falsePos[ii][jj]-falseNeg[ii][jj]);
	mmissCut = ii;
	pmissCut = jj;
      }
      
    }//end jj
  }//end ii
  cout<<"optimal pmiss cut: "<<pmiss_cut[pmissCut]<<" optimal mmiss cut: "<<mmiss_val[mmissCut]<<" false Pos: "<<falsePos[mmissCut][pmissCut]<<" false Neg: "<<falseNeg[mmissCut][pmissCut]<<endl;

  
  TGraph *grP[7];
  TGraph *grN[7];
  TMultiGraph *mgP = new TMultiGraph();
  TMultiGraph *mgN = new TMultiGraph();
  TLegend* legendP = new TLegend(0.7,0.6,0.8,0.8);
  TLegend* legendN = new TLegend(0.7,0.6,0.8,0.8);
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

  //f->Close();
  double totEvents;
  canvas->Update();

  mgP->Draw("a");
  legendP->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  mgN->Draw("a");
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

  h_emiss->SetLineColor(kRed);
  h_emissSmear->SetLineColor(kBlue);
  h_emiss->Draw();
  h_emissSmear->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom->SetLineColor(kRed);
  h_nmomSmear->SetLineColor(kBlue);
  h_nmom->Draw();
  h_nmomSmear->Draw("same");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_q2->Draw();
  totEvents = h_q2->GetEntries();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_smear->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  hist_pq_thetapq->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  hist_pmiss_xb->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpmSm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss->SetLineColor(kRed);
  h_mmissSmear->SetLineColor(kBlue);
  h_mmiss->Draw();
  h_mmissSmear->Draw("same");
  //h_zvtx_precut->Scale(1./h_zvtx_precut->GetMaximum());
  //h_zvtx_cut->Scale(1./h_zvtx_cut->GetMaximum());
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_zvtx_precut->Draw();
  canvas->Print( (pdf_file_name + ")").c_str());





  fout->Write();
  fout->Close();
  cout<<"Good events: "<<totEvents<<endl;

}
