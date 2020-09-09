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

//Define some globale functions and constants
double accMin = 0.4;
double sq(double x){
  return x*x;
}

void help_message()
{
  cerr<< "Argumets: ./get_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Nucleus A for Ratio] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-v: Verbose\n\n";
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
  //Get target info
  int A = atoi(argv[3]);
  target_Info targInfo(A);  
  target_Info secondTargInfo(atoi(argv[4]));
  bool verbose = false;
  
  int c;
  while ((c=getopt (argc-4, &argv[4], "hv")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
      case 'v':
	verbose = true;
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cout<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");
  double binxB[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.58,1.77,2};
  int numbinxB = ( sizeof(binxB)/sizeof(binxB[0]) ) - 1;

  //Make a histogram list to make things easier
  vector<TH1*> hist_list;
  TH1D * hist_xB =  new TH1D("hist_xB" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",20,0,5);
  hist_list.push_back(hist_QSq);
  TH1D * hist_pp = new TH1D("hist_pp","hist;pp;Counts",40,-0.6,0.6);
  hist_list.push_back(hist_pp);
  TH1D * hist_pt = new TH1D("hist_pt","hist;pt;Counts",20,0,0.6);
  hist_list.push_back(hist_pt);

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
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);
    
    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    

    TVector3 ve = myInfo.getVector(0);
    double ThetaE = ve.Theta() * 180 / M_PI;
    TVector3 vq = vBeam - ve;
    double omega = vBeam.Mag() - ve.Mag();
    double phi = ve.Phi() * 180 / M_PI;
    TVector3 vLead = myInfo.getVector(1);
    TVector3 vMiss = vLead - vq;
    double mMiss = myInfo.getMassMiss(1);
    double poq = myInfo.getPoQ(1);
    double thetapq = myInfo.getThetaPQ(1);
    double thetapMq = vq.Angle(vMiss) * (180/M_PI);
    double pMiss_p = vMiss.Dot(vq.Unit());
    double  pMiss_t = vMiss.Perp(vq);

    double weight = 1;
    
    //Only take particles that pass fiducials and acceptances for both targets
    if(targInfo.e_acc(ve) < accMin){continue;}
    if(targInfo.p_acc(vLead) < accMin){continue;}
    if(!targInfo.pass_semi_fid(ve,vLead)){ continue; }      
    if(secondTargInfo.e_acc(ve) < accMin){continue;}
    if(secondTargInfo.p_acc(vLead) < accMin){continue;}
    if(!secondTargInfo.pass_semi_fid(ve,vLead)){ continue; }

    //Establish vertex Cuts
    if(!targInfo.semiVTXInRange(myInfo,1)){
      continue;
    }


    //Apply the corrections to the cross sections
    // Acceptance correction
    // Transparancy correcton
    // Radiative correction
    // Normalize with luminosity and A
    weight = weight / targInfo.semi_acc(ve,vLead);
    weight = weight / targInfo.getTrans();
    weight = weight / targInfo.getRadCorr(ThetaE,xB);
    weight = weight / (targInfo.getLum() * A);


    hist_xB->Fill(xB,weight);	  	        
    hist_QSq->Fill(QSq,weight);
    hist_pp->Fill(pMiss_p,weight);
    hist_pt->Fill(pMiss_t,weight);

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
