#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring> 

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"

#include "constants.h"

using namespace std;

const TVector3 vBeam(0.,0.,4.461);
//Define some useful functions
double sq(double x){ return x*x; };

int getSec(const TVector3 ve);

int main(int argc, char ** argv){

  if( argc != 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "nMass /path/to/input/data/file /path/to/output/root/file \n";

    return -1;
    
  }
  
  //Get A trees. Open output file.
  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");

  //Make trees and histograms for data
  TTree * inTree = (TTree*)inputFile->Get("T");

  //Define constants that I should redifine elsewhere
  const int pCode = 2212;
  const int nCode = 2112;
  const double m_3He = 3.016 * mU;

  cerr<<"Nucleus file has been opened from input file. \n";

  
  //Define varialbes needed for input tree and make tree
  Int_t nPar;
  Int_t parID[19];
  Double_t px[19],py[19],pz[19];
  inTree->SetBranchAddress("nParticles",&nPar);
  inTree->SetBranchAddress("Part_type",parID);
  inTree->SetBranchAddress("mom_x",px);
  inTree->SetBranchAddress("mom_y",py);
  inTree->SetBranchAddress("mom_z",pz);

  //Make a vector of histograms to loop over them
  vector<TH1*> hist_list;
  
  // Make histogram for all sectors
  TH1D * hist_nMass = new TH1D("nMass_diff","hist_nMass_diff;nMass_diff;Counts",100,-0.5,0.5);
  hist_list.push_back(hist_nMass);
  
  //Make histogram for each sector
  TH1D * hist_nMass_bySec[6];
  for(int i=0; i<6; i++){
    char temp[100];
    sprintf(temp,"nMass_diff_%d",i);
    hist_nMass_bySec[i] = new TH1D(temp,"nMass_diff_by_Sec;nMass_diff;Counts",100,-0.5,0.5);
    hist_list.push_back(hist_nMass_bySec[i]);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();  
  }
  

  cerr<<"input TTree opened from input file and output histogram created. \n";

 int fin = inTree->GetEntries();
  cout<<fin<<"\n";
  
  //Loop over TTree
  for(int i = 0; i < fin; i++){

    inTree->GetEntry(i);
    
    //Get Momentum
    TVector3 ve(px[0],py[0],pz[0]);
    TVector3 vp1(px[1],py[1],pz[1]);
    TVector3 vp2(px[2],py[2],pz[2]);
    TVector3 vn = vBeam - (ve + vp1 + vp2);

    //Get Energy
    double Ep1 = sqrt(vp1.Mag2()+sq(mP));
    double Ep2 = sqrt(vp2.Mag2()+sq(mP));
    double En = (m_3He + vBeam.Mag()) - (ve.Mag() + Ep1 +Ep2);


    //if((sq(En) - vn.Mag2()) < 0) continue;
    
    double nMass = sqrt(sq(En) - vn.Mag2());
    double nMass_diff=nMass-mN;

    double theta = ve.Theta() * 180 /M_PI;
    if(theta > 16) continue;
    
    //Fill Histogram
    hist_nMass->Fill(nMass_diff,1);
    hist_nMass_bySec[getSec(ve)]->Fill(nMass_diff,1);
  }


  
  //Now write out the histograms and close files
  outputFile->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();  
  }
  outputFile->Close();
  inputFile->cd();
  inputFile->Close();
  

  return 0;
}

int getSec(const TVector3 ve){

  double phi = ve.Phi() * 180 / M_PI;

  if(phi<=-150) return 0;
  else if(phi<=-90) return 1;
  else if(phi<=-30) return 2;
  else if(phi<=30) return 3;
  else if(phi<=90) return 4;
  else if(phi<=150) return 5;
  else return 0; 

}
