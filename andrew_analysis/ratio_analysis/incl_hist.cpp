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
  cerr<< "Argumets: ./incl_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Nucleus A for Ratio] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-n: Set a custom minimum to the Z vertex [cm]\n"
      <<"-x: Set a custom maximum to the Z vertex [cm]\n\n";
}


int main(int argc, char ** argv){
     
  if (argc < 2)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    } 
  if (strcmp(argv[1], "-h")==0)
    {
      help_message();
      return -1;
    }
  if (argc < 5)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    }

  //Read in arguments

   TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");
  int A = atoi(argv[3]);
  target_Info secondTargInfo(atoi(argv[4]));
  //Get target info
  target_Info targInfo(A);  
  bool verbose = true;
  
  int c;
  while ((c=getopt (argc-4, &argv[4], "hn:x:")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
      case 'n':
	targInfo.change_vtxMin(atof(optarg));
	break;
      case 'x':
	targInfo.change_vtxMax(atof(optarg));
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cout<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");
  double finebinxB[] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.40,1.50,1.62,1.76,1.88,2};
  int finenumbinxB =  ( sizeof(finebinxB)/sizeof(finebinxB[0]) ) - 1;

  //Make a histogram list to make things easier
  vector<TH1*> hist_list;
  TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_xB_noRad =  new TH1D("hist_xB_noRad_Incl" ,"hist;xB;Counts",finenumbinxB,finebinxB);
  hist_list.push_back(hist_xB_noRad);
  TH1D * hist_xB =  new TH1D("hist_xB_Incl" ,"hist;xB;Counts",finenumbinxB,finebinxB);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",40,0,5);
  hist_list.push_back(hist_QSq);

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
  }

  
  cerr<<"Histograms and Trees successfully created\n";
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Set addresses for inTree
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

  int counter = 0;
  //Loop over TTree
  for(int i = 0; i < inTree->GetEntries(); i++){
    double weight = 1;
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);

    TVector3 ve = myInfo.getVector(0);
    TVector3 vq = vBeam - ve;
    double omega = vBeam.Mag() - ve.Mag();
    double phi = ve.Phi() * 180 / M_PI;
    double theta = ve.Theta() * 180 / M_PI;
    
    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    
    if(!targInfo.evtxInRange(vtxZCorr[0])){continue;}
    if(targInfo.e_acc(ve) < accMin){continue;}
    if(secondTargInfo.e_acc(ve) < accMin){continue;}
    if(QSq < 1){continue;}
    if(xB < 0.15){continue;}
    if(xB > 2.001){continue;}
    double eff = 1;
    eff = targInfo.incl_acc(ve);
    weight = weight / eff;
    weight = weight/(targInfo.getLum() * A);

    hist_xB_noRad->Fill(xB,weight);	  	  

    weight = weight / targInfo.getRadCorr(theta,xB);
    hist_xB->Fill(xB,weight);	  	      
    hist_Cross->Fill(1,1);
    hist_QSq->Fill(QSq,weight);
    
    counter++;
    	
  }

  cout<<"There were "<<counter<<" events written\n";
  cout<<"The total luminosity for these runs is: "<<targInfo.getLum()<<"\n";
  cerr<<"Finished filling histogram\n";
   
  inputFile->cd();
  inputFile->Close();

  
  //Now write out
  outputFile->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  outputFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
