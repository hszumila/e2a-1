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

const double Me      = .000511;         // Electron     mass in GeV
const double Mp      = 0.93827;         // Proton       mass in GeV
const double Mn      = 0.93957;         // Neutron      mass in GeV
const double He3_be  = 0.00770;         // 3He binding energy in GeV
const double M_tar   = 2*Mp+Mn-He3_be;  // 3Helium  mass in GeV (target)

int main(){

  gROOT->SetBatch(true);
  //gStyle->SetOptStat(0);

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
  TFile *f = new TFile("../../../data/n_eff/skim_He3_pp.root");
  TFile *fout = new TFile("output/neff_He3.root","RECREATE");
  
  TCanvas *canvas = new TCanvas("canvas","output plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/neff_He3.pdf";

  //load branches
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

  //loop events
  for (int event =0 ; event < intree->GetEntries() ; event++)
    {
      if (event%10000==0)
      	cout << "Working on event " << event <<endl;
      
      intree->GetEvent(event);

      // Get the electron and q vectors
      TVector3 mom_e0(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e0.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e0;

      //plot the EC time resolution
      h_tres->Fill(sc_time[0]-ec_time[0]-(sc_path[0]-ec_path[0])/30.);

      // Loop over particles looking for exclusive 3He(e,e'ppn)
      int leadingID = -1;
      if (nParticles != 4) {continue;}

      int p1 = -1;
      int p2 = -1;
      int n0 = -1;
      for (int part=1 ; part < nParticles ; part++)
	{
	  //check if vertex in 3He
      	  //if (vtxZ[part]<-2.5 || vtxZ[part]>-0.5){continue;}
	  
	  if (type[part] == 2112){n0 = part;}
	  else if (type[part] == 2212){
	    if (p1==-1){p1=part;}
	    else{p2=part;}
	  }

	}

      if (p1<0 || p2<0 || n0<0){continue;}

      TLorentzVector mom_e,mom_p1,mom_p2, mom_n, E0, V4_tar;
      TVector3 P_miss_eppn;

      E0.SetXYZT(0.,0.,sqrt(Ebeam*Ebeam-Me*Me),Ebeam);
      V4_tar.SetXYZT(0.,0.,0.,M_tar);

      mom_e.SetXYZM(mom_x[0],mom_y[0],mom_z[0],Me);// Electron
      mom_n.SetXYZM(mom_x[n0],mom_y[n0],mom_z[n0],Mn);// Electron
      mom_p1.SetXYZM(mom_x[p1],mom_y[p1],mom_z[p1],Mp);// Proton 1
      mom_p2.SetXYZM(mom_x[p2],mom_y[p2],mom_z[p2],Mp);// Proton 2

      // Calculate missing mass and momentum	
      double M_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2).M();
      P_miss_eppn = (E0 + V4_tar - mom_e - mom_p1 - mom_p2 ).Vect();

      if (M_miss_eppn<0.8 || M_miss_eppn>1.05){continue;}
      
      //calc pmiss
      double angle = P_miss_eppn.Angle(mom_n.Vect())*180.0/M_PI;
      if (angle>10.){continue;}
      
      h_rmDot->Fill(angle);
      h_nmom->Fill(mom_n.P());
      h_pmom->Fill(mom_p1.P());
      h_pmom->Fill(mom_p2.P());
      h_miss->Fill(P_miss_eppn.Mag());
      h_mmiss->Fill(M_miss_eppn);

      double pmiss_diff = mom_n.P()-P_miss_eppn.Mag();
      h_pdiff_tot->Fill(pmiss_diff);
      h_pdiffVp->Fill(pmiss_diff,mom_n.P());
      h_pmVpm->Fill(mom_n.P(),P_miss_eppn.Mag());

      
      //plot the p_measured - p_miss in slices
      if(mom_n.P()>=0.5 && mom_n.P()<1){
	h_pdiff_a->Fill(pmiss_diff);
	nmom_a->Fill(mom_n.P());
      }
      if (mom_n.P()>=1 && mom_n.P()<1.4){
	h_pdiff_b->Fill(pmiss_diff);
	nmom_b->Fill(mom_n.P());
      }
      if (mom_n.P()>=1.4 && mom_n.P()<2.){
	h_pdiff_c->Fill(pmiss_diff);
	nmom_c->Fill(mom_n.P());
      }
      
    }//end loop

  
  canvas->Update();


  h_tres->Draw();
  TF1 *g1 = new TF1("g1","gaus",-0.5,1);
  //2261:
  //TF1 *g1 = new TF1("g1","gaus",1,2);
  h_tres->Fit(g1,"R");
  double delt = abs(g1->GetParameter(2));
  //g1->Draw("same")
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pdiff_tot->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pdiff_a->Draw();
  TF1 *g2 = new TF1("g2","gaus",-0.5,0.5);
  h_pdiff_a->Fit(g2,"R");
  double sigA = abs(g2->GetParameter(2));
  double statErrA = 1./sqrt(h_pdiff_a->Integral());
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pdiff_b->Draw();
  TF1 *g3 = new TF1("g3","gaus",-0.2,0.5);
  h_pdiff_b->Fit(g3,"R");
  double sigB = abs(g3->GetParameter(2));
  double statErrB = 1./sqrt(h_pdiff_b->Integral());
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pdiff_c->Draw();
  TF1 *g4 = new TF1("g4","gaus",-0.2,0.5);
  //2261:
  //TF1 *g4 = new TF1("g4","gaus",0,1);
  h_pdiff_c->Fit(g4,"R");
  double sigC = abs(g4->GetParameter(2));
  double statErrC = 1./sqrt(h_pdiff_c->Integral());
  canvas->Print( (pdf_file_name + "(").c_str());


  nmom_a->Draw();
  double meanA = nmom_a->GetMean();
  double rmsA = nmom_a->GetRMS();
  canvas->Print( (pdf_file_name + "(").c_str());

  nmom_b->Draw();
  double meanB = nmom_b->GetMean();
  double rmsB = nmom_b->GetRMS();
  canvas->Print( (pdf_file_name + "(").c_str());

  nmom_c->Draw();
  double meanC = nmom_c->GetMean();
  double rmsC = nmom_c->GetRMS();
  canvas->Print( (pdf_file_name + "(").c_str());
 
  h_pdiffVp->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pmVpm->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
	
  h_rmDot->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_pmom->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_miss->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  TMultiGraph *mg = new TMultiGraph();

  //double mom[3] = {meanA,meanB,meanC};
  //double mombin[3] = {rmsA,rmsB,rmsC};
  double mom[3] = {0.75,1.2,1.7};
  double mombin[3] = {0.25,0.2,0.3};
  double momUncert[3] = {sigA/meanA, sigB/meanB, sigC/meanC};
  double momUncertErr[3] = {sigA*statErrA, sigB*statErrB, sigC*statErrC};//statistical
  
  TGraphErrors *tg =new TGraphErrors(3,mom,momUncert,mombin,momUncertErr);
  tg->SetMarkerColor(kBlue);
  tg->SetMarkerStyle(22);
  tg->SetLineColor(kBlue);
  tg->SetLineWidth(2);
  mg->Add(tg);
  mg->Draw("ap");

 
  
  double xcal = 516.9152; //ec dist in cm
  double Mn = 0.939565413;
  double pres[10];
  double pp[10];
  for (int ii=0; ii<10; ii++){
    pp[ii] = (double) ii/4.0;
    pres[ii] = pp[ii]/pow(Mn,2.0)*sqrt(pow(Mn,2.0)+pow(pp[ii],2.0))*delt/xcal*30.0;
  }

  TF1 *ffit = new TF1("ffit",Form("sqrt(pow(x/pow(%f,2.0)*sqrt(pow(%f,2.0)+pow(x,2.0))*%f/%f*30.0,2.0)+pow([0],2.0))",Mn,Mn,delt,xcal),0,2.1);
  ffit->SetParameter(0,1);
  mg->Fit(ffit,"R");
  ffit->SetLineColor(kBlue);
  ffit->Draw("lsame");
  
  mg->GetYaxis()->SetTitle("#Delta p/p");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->SetRangeUser(0.0,3);
  mg->GetYaxis()->SetRangeUser(0.0,0.3);
  mg->GetXaxis()->SetTitle("Neutron momentum resolution");
  gPad->Modified();

  TGraph *g_est = new TGraph(10,pp,pres);

  g_est->SetLineColor(kBlack);
  g_est->SetLineWidth(2.5);
  g_est->Draw("lsame");
  canvas->Print( (pdf_file_name + ")").c_str());

  
}
