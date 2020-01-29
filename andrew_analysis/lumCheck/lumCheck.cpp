#include <fstream>
#include <iostream>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

#include "target_Info.h"

using namespace std;


int main(int argc, char ** argv){

  if( argc != 5){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "lumCheck /path/to/text/file /path/to/root/data/files[path/incl_skim_%i.root] /path/to/output/root/file [A]\n";

    return -1;
    
  }

  ifstream runlist;
  runlist.open(argv[1]);
  char * filepath = argv[2];
  TFile * outputFile = new TFile(argv[3],"RECREATE");
  int A = atoi(argv[4]);
  //Get target info for the vertext cuts
  target_Info myInfo(A);
  
  //Get passed first six lines
  char trash[256];  
  for (int i=0;i<6;i++){
  runlist .getline (trash,256);
  }
  
  //Define variables to pull from text file
  int runNum;
  double charge;
  double livecharge;
  double triggers;

  //Make vectors to put into a TGraph
  TGraph * forwardGraph = new TGraph();
  forwardGraph->SetName("luminosityForward");
  TGraph * backwardGraph = new TGraph();
  backwardGraph->SetName("luminosityBackward");

  while( (!runlist.eof()) && (runlist.is_open()==true) ){

    //Get the run number and livecharge from the text file for each run
    runlist >> runNum >> charge >> livecharge >> triggers;

    //Open up the file
    TFile* f1 = new TFile(Form(filepath,runNum));
    if(f1->IsZombie()){
      cout<<"Could not open file: "<<Form(filepath,runNum)<<"\n Skipping file\n";
      continue;
    }
    TTree* t1 = (TTree*)f1->Get("T");

    //Make histograms for the forward and backward kinematics
    TH1D * h1_xB_for = new TH1D(Form("xB_%i_for",runNum),"xB_for;xB;counts",50,1.,2.);
    TH1D * h1_xB_back = new TH1D(Form("xB_%i_back",runNum),"xB_back;xB;counts",50,1.,2.);
    
    //Number of events in forward and backward direction
    double numEventsFor = 0;
    double numEventsBack = 0;
    
    //Define variables for tree that I will need to cut on
    double Xb,Q2,vtxCorr;
    double px[5],py[5],pz[5];
    t1->SetBranchAddress("Xb",&Xb);
    t1->SetBranchAddress("Q2",&Q2);
    t1->SetBranchAddress("mom_x",px);
    t1->SetBranchAddress("mom_y",py);
    t1->SetBranchAddress("mom_z",pz);
    t1->SetBranchAddress("vtx_z_cor",&vtxCorr);
    //Loop over TTree
    for(int i = 0; i < t1->GetEntries(); i++){
      t1->GetEntry(i);
      


      //Get varibles needed for cuts
      TVector3 ve(px[0],py[0],pz[0]);
      double theta = ve.Theta() * 180 / M_PI;
      double p = ve.Mag();

      //if(!myInfo.evtxInRange(vtxCorr)){
      //	continue;
      //}
      
      //Check histograms and values in the forward and backward direction
      if((theta>18) && (theta<22)){
	if((p>2.2) && (p<3)){
	  numEventsFor = numEventsFor + 1;
	  h1_xB_for->Fill(Xb);
	}
      }
      else if((theta>22) && (theta<28)){
	if((p>1.5) && (p<2.2)){
	  numEventsBack = numEventsBack + 1;
	  h1_xB_back->Fill(Xb);
	}
      }   
    }
    
    //write out histograms
    outputFile->cd();
    h1_xB_for->Scale(1/livecharge);
    h1_xB_for->SetMaximum(0.003);
    h1_xB_for->Write();
    h1_xB_back->Scale(1/livecharge);
    h1_xB_back->SetMaximum(0.0002);
    h1_xB_back->Write();

    //Make vectors of results for the number of events
    forwardGraph->SetPoint(forwardGraph->GetN(),runNum,(numEventsFor/livecharge));
    backwardGraph->SetPoint(backwardGraph->GetN(),runNum,(numEventsBack/livecharge));
       
    //Now close out the file for this runNumber
    f1->cd();
    f1->Close();
    
    cout << Form(filepath,runNum) <<"\n Forward: "<<(numEventsFor/livecharge)<<"      Backward: "<<(numEventsBack/livecharge)<<endl;
     
  }
  
  forwardGraph->SetMarkerStyle(8);
  backwardGraph->SetMarkerStyle(8);

  //Make the Tgraphs and write out
  outputFile->cd();
  forwardGraph->Write();
  backwardGraph->Write();
  outputFile->Close();
  
return 0;
}

