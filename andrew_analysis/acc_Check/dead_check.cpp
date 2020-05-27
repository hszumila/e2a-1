#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 
#include <unistd.h> 
#include <iterator>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"

#include "target_Info.h"
#include "event_Info.h"

using namespace std;
double accMin = 0.4;
double sq(double x){
  return x*x;
}

void help_message()
{
  cerr<< "Argumets: ./incl_hist /path/to/output/file [Nucleus A] /path/to/input/skim/files... \n\n";
}


int main(int argc, char ** argv){
     
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    } 

  //Read in arguments
  target_Info targInfo(atoi(argv[2]));  
  char temp[100];

  TH2D * hist_1e[argc-3];
  TH2D * hist_2e[argc-3];
  TH2D * hist_4e[argc-3];
  for(int i = 0; i < (argc-3); i++){
    sprintf(temp,"Hist_%d_energy_%d",i,1500);
    hist_1e[i] = new TH2D(temp,"Hist;phi;theta",360,-180,180,50,10,60);

    sprintf(temp,"Hist_%d_energy_%d",i,2000);
    hist_2e[i] = new TH2D(temp,"Hist;phi;theta",360,-180,180,50,10,60);

    sprintf(temp,"Hist_%d_energy_%d",i,2500);
    hist_4e[i] = new TH2D(temp,"Hist;phi;theta",360,-180,180,50,10,60);
  }

  for(int i = 3; i < argc; i++){
    TFile * inputFile = new TFile(argv[i]);
    cout<<argv[i]<<" has been opened\n";
    //Set addresses for inTree
    TTree * inTree = (TTree*)inputFile->Get("T");
    //Define variables needed for histograms and tree
    Int_t nPar;
    Int_t parID[19];
    Double_t xB,QSq;
    Double_t momx[19],momy[19],momz[19],vtxZCorr[19];
    const Double_t Ebeam = 4.461;
    inTree->SetBranchAddress("nParticles",&nPar);
    inTree->SetBranchAddress("Xb",&xB);
    inTree->SetBranchAddress("Q2",&QSq);  
    inTree->SetBranchAddress("Part_type",parID);
    inTree->SetBranchAddress("mom_x",momx);
    inTree->SetBranchAddress("mom_y",momy);
    inTree->SetBranchAddress("mom_z",momz);
    inTree->SetBranchAddress("vtx_z_cor",vtxZCorr);

    for(int j = 0; j < inTree->GetEntries(); j++){
      inTree->GetEntry(j);
      event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);
      TVector3 ve = myInfo.getVector(0);
      double eE = ve.Mag();
      TVector3 vq = vBeam - ve;
      double omega = vBeam.Mag() - ve.Mag();
      double phi = ve.Phi() * 180 / M_PI;
      phi+=30;
      if(phi>180){phi-=360;}
      double theta = ve.Theta() * 180 / M_PI;

      //Do some things for all events
      if((j%5000) == 0){cout << (j*100.)/(inTree->GetEntries()) <<"% complete \n";}
    
      if(nPar > 19){
	cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
	return -1;
      } 

      //Fiducial Cuts
      if(targInfo.e_acc(ve) < accMin){continue;}
      if(!targInfo.pass_incl_fid(ve)){continue;}

      if(xB < 0.15){continue;}
      if(xB > 2){continue;}

      double weight = 1;
      //weight = weight / targInfo.incl_acc(ve);
      //weight = weight / targInfo.getRadCorr(theta,xB);

      if((xB > 0.6) && (xB < 0.7)){
	hist_1e[i-3]->Fill(phi,theta,weight);      
      }
      if((xB > 0.9) && (xB < 1)){
	hist_2e[i-3]->Fill(phi,theta,weight);      
      }
      if((xB > 1.2) && (xB < 1.3)){
	hist_4e[i-3]->Fill(phi,theta,weight);      
      }

    }
   
    inputFile->cd();
    inputFile->Close();
  }  
  cerr<<"Finished filling histogram\n";

  //Now write out
  TFile * outputFile = new TFile(argv[1],"RECREATE");
  for(int k=0; k<(argc-3); k++){
    cout<<k<<"\n";
    outputFile->cd();
    hist_1e[k]->Write();
    hist_2e[k]->Write();
    hist_4e[k]->Write();
  }
  outputFile->Close();
  cerr<< argv[1]<<" has been completed. \n\n\n";
  
  return 0;
}
