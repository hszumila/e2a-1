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

int getGroup(int spot, int total){
  double x = (double)spot/(double)total;
  if((x>=0) && (x<0.25)){return 0;}
  else if((x>=0.25) && (x<0.5)){return 1;}
  else if((x>=0.5) && (x<0.75)){return 2;}
  else if((x>=0.75) && (x<1)){return 3;}
  return 3;
}
void help_message()
{
  cerr<< "Argumets: ./get_hist /path/to/output/file [Nucleus A] /path/to/input/skim/files...  \n\n\n\n";
}


int main(int argc, char ** argv){
     
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    } 

  TFile * outputFile = new TFile(argv[1],"RECREATE");
  int A = atoi(argv[2]);
  target_Info targInfo(A);  
  int counter = 0;

  double binT[] = {14,14.25,14.5,14.75,15,15.25,15.5,15.75,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,26,27,28,29,30,32,34,36,38,40};
  int numbinT = ( sizeof(binT)/sizeof(binT[0]) ) - 1;

  TH1D * hist_T1[4][6];
  TH2D * hist_T2[4][6];
  char temp[100];
  for(int m = 0; m<1; m++){
    for(int i = 0; i<6; i++){
      sprintf(temp,"hist_Theta_grp%d_sec%d",m,i);
      hist_T1[m][i] = new TH1D(temp,"Theta Distribution;Theta;counts",numbinT,binT);
      hist_T1[m][i]->Sumw2();
      hist_T2[m][i] = new TH2D(temp,"Theta Distribution;Theta;xB;counts",numbinT,binT,20,0.8,2);
      hist_T2[m][i]->Sumw2();
    }
  }


  TH1D * hist_sec[argc-3];
  for(int i = 0; i<(argc-3); i++){
    sprintf(temp,"hist_sec_run%d",i);
    hist_sec[i] = new TH1D(temp,"sec distribution; sector; counts",6,-0.5,5.5);
    hist_sec[i]->Sumw2();
  }
  
  cerr<<"Histograms and Trees successfully created\n";

  for(int k = 3; k < argc; k++){
    TFile * inputFile = new TFile(argv[k]);
    cout<<"File "<< argv[k] <<" has been opened\n";

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Set addresses for inTree
    //Define variables needed for histograms and tree
    TTree * inTree = (TTree*)inputFile->Get("T");
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
    int myGroup = getGroup(k-3,argc-3);    
    //Loop over TTree
    for(int i = 0; i < inTree->GetEntries(); i++){
      inTree->GetEntry(i);
      event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);     
      TVector3 ve = myInfo.getVector(0);
      double phi = ve.Phi() * 180 / M_PI;
      double theta = ve.Theta() * 180 / M_PI;
      int eventSector = targInfo.getSecPhi(phi);
      //Do some things for all events
      if((i%5000) == 0){
	cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
      }

      if(targInfo.e_acc(ve) < accMin){continue;}
      if(!targInfo.pass_incl_fid(ve)){continue;}      
      if(!targInfo.eVTXInRange(vtxZCorr[0])){continue;}
      
      if((xB<1.4) || (xB>1.8)){continue;}

      double weight = 1;
      counter++;
      hist_sec[k-3]->Fill(eventSector,weight);
      //myGroup
      hist_T1[0][eventSector]->Fill(theta,weight);
      hist_T2[0][eventSector]->Fill(theta,xB,weight);
    }

    inputFile->cd();
    inputFile->Close();
  }
  cout<<"There were "<<counter<<" events written\n";

  //Now write out
  outputFile->cd();
  for(int m = 0; m<1; m++){
    for(int j = 0; j < 6; j++){
      hist_T1[m][j]->Write();
      hist_T2[m][j]->Write();
    }   
  }
  for(int i = 0; i<(argc-3); i++){
    hist_sec[i]->Write();
  }
  outputFile->Close();
  cerr<< argv[1]<<" has been completed. \n\n\n";
  
  return 0;
}
