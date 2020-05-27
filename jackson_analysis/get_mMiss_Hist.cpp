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
  cerr<< "Argumets: ./get_mMiss_Hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-s: Input a sector\n\n";
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
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    }

  //Read in arguments

  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");
  int A = atoi(argv[3]);
  target_Info targInfo(A);  
  bool verbose = true;
  bool doOneSec = false;
  int mySec = -1;
  int c;
  while ((c=getopt (argc-4, &argv[4], "hs:")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
      case 's':
	doOneSec=true;
	mySec=atof(optarg);
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cout<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");

  TH1D * h_mMiss = new TH1D("ep_mMiss","ep;mMiss [GeV];Counts",30,0.5,1.1);
  h_mMiss->Sumw2();
  TH2D * h_mMiss_Pm = new TH2D("ep_mMiss_Pm","ep;mMiss [GeV]; pmiss [GeV];Counts",30,0.5,1.1,24,0.4,1.0);
  h_mMiss_Pm->Sumw2();
  TH2D * h_mMiss_QSq = new TH2D("ep_mMiss_QSq","ep;mMiss [GeV]; QSq [GeV^2];Counts",30,0.5,1.1,32,1.0,5.0);
  h_mMiss_QSq->Sumw2();
  TH2D * h_mMiss_xB = new TH2D("ep_mMiss_xB","ep;mMiss [GeV]; xB;Counts",30,0.5,1.1,32,1.2,2.0);
  h_mMiss_xB->Sumw2();

  
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
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);


    double mMiss = myInfo.getMassMiss(1);
    TVector3 ve = myInfo.getVector(0);
    TVector3 vq = vBeam - ve;
    TVector3 vLead = myInfo.getVector(1);
    TVector3 vMiss = vq - vLead;
    double mA = targInfo.getMass();

    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    

    
    //Only take particles that pass fiducials and acceptances
    if(targInfo.e_acc(ve) < accMin){continue;}
    if(targInfo.p_acc(vLead) < accMin){continue;}
    if(!targInfo.pass_semi_fid(ve,vLead)){ continue; }
      
    //Establish vertex Cuts
    //You can look at these histograms if you want to see how your vertex cuts look
    //hist_eVTX->Fill(vtxZCorr[0],1);
    //hist_pVTX->Fill(vtxZCorr[1],1);
    //hist_vtx->Fill(vtxZCorr[0],vtxZCorr[1],1);
    if(!targInfo.semiFixedVTXInRange(myInfo,1)){
      continue;
    }


    //Apply the corrections to the cross sections
    // Acceptance correction
    // Transparancy correcton
    // Radiative correction
    // Normalize with luminosity and A
    double weight = 1;
    //weight = weight / targInfo.semi_acc(ve,vLead);
    //weight = weight / targInfo.getTrans();
    //weight = weight / targInfo.getRadCorr(theta,xB);
    //weight = weight / (targInfo.getLum() * A);

    counter++;

    h_mMiss->Fill(mMiss,weight);
    h_mMiss_Pm->Fill(mMiss,vMiss.Mag(),weight);
    h_mMiss_QSq->Fill(mMiss,QSq,weight);
    h_mMiss_xB->Fill(mMiss,xB,weight);    

  }

  cout<<"There were "<<counter<<" events written\n";
  cout<<"The total luminosity for these runs is: "<<targInfo.getLum()<<"\n";
  cerr<<"Finished filling histogram\n";
   
  inputFile->cd();
  inputFile->Close();

  
  //Now write out
  outputFile->cd();
  h_mMiss->Write();
  h_mMiss_Pm->Write();
  h_mMiss_QSq->Write();
  h_mMiss_xB->Write();
  outputFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
