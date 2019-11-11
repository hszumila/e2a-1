#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <math.h> 
#include <unistd.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"

#include "target_Info.h"
#include "constants.h"
#include "Acceptance.h"

using namespace std;

//Define some globale functions and constants
double sq(double x){
  return x*x;
}
int getSec(const double phi);

void help_message()
{
  cerr<< "Argumets: ./incl_hist /path/to/input/skim/file [Nucleus a] [minimum momentum in GeV] [maximum momentum in GeV] [minimum theta in degrees] [maximum theta in degrees] [minimum phi in degrees] [maximum phi in degrees] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-v: Verbose\n"
      <<"-n: Set a custom minimum to the Z vertex [cm]\n"
      <<"-x: Set a custom maximum to the Z vertex [cm]\n"
      <<"-a: Keep only points will low acceptance\n\n";
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
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    }

  //Read in arguments
  TFile * inputFile = new TFile(argv[1]);
  int A = atoi(argv[2]);
  double momMin = atoi(argv[3]);
  double momMax = atoi(argv[4]);
  double theMin = atoi(argv[5]);
  double theMax = atoi(argv[6]);
  double phiMin = atoi(argv[7]);
  double phiMax = atoi(argv[8]);
  //Get target info
  target_Info targInfo(A);
  
  bool verbose = false;
  bool lowAcc = false;
  
  int c;
  while ((c=getopt (argc-8, &argv[8], "hvn:x:a")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
      case 'v':
	verbose = true;
	break;
      case 'n':
	targInfo.change_vtxMin(atof(optarg));
	break;
      case 'x':
	targInfo.change_vtxMax(atof(optarg));
	break;
      case 'a':
	lowAcc = true;
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cerr<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");
  //Set addresses for inTree
  //Define variables needed for histograms and tree
  Int_t nPar;
  Int_t parID[19];
  Double_t xB,QSq;
  Double_t px[19],py[19],pz[19],vtxZCorr[19];
  const Double_t Ebeam = 4.461;
  inTree->SetBranchAddress("nParticles",&nPar);
  inTree->SetBranchAddress("Xb",&xB);
  inTree->SetBranchAddress("Q2",&QSq);  
  inTree->SetBranchAddress("Part_type",parID);
  inTree->SetBranchAddress("mom_x",px);
  inTree->SetBranchAddress("mom_y",py);
  inTree->SetBranchAddress("mom_z",pz);
  inTree->SetBranchAddress("vtx_z_cor",vtxZCorr);

  //Loop over TTree
  for(int i = 0; i < inTree->GetEntries(); i++){

    inTree->GetEntry(i);    
    //Display completed
    if(((i%100000) == 0) && verbose){
      cerr << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    if(nPar > 19){
      cout<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    
    //Only look in a range of QSq
    if((QSq<1) || (QSq>=5)){ 
      continue;
    }
    //First get rid of events outside of the vertex
    if(!targInfo.evtxInRange(vtxZCorr[0])){
      continue;
    }

    TVector3 ve(px[0],py[0],pz[0]);
    double pe = ve.Mag();
    double theta = ve.Theta() * 180 / M_PI;
    double phi = ve.Phi() * 180 / M_PI;
    int sec = getSec(phi);

    //Only look in at data with a specific momentum value
    if( (pe < momMin ) || (pe > momMax) ){
      continue;
    } 
    if( (theta < theMin ) || (theta > theMax) ){
      continue;
    } 
    if( (phi < phiMin ) || (phi > phiMax) ){
      continue;
    } 

    //Only write out one if the acceptance is low
    if(lowAcc && (targInfo.e_acc(ve)>0.1)){
      continue;
    }
    
    cout<<pe<<" "<<theta<<" "<<phi<<" "<<xB<<" "<<QSq<<" "<<targInfo.e_acc(ve)<<" \n";
  
  }
  
  inputFile->cd();
  inputFile->Close();
  return 0;
}

int getSec(const double phi){
  

  double nphi=phi;//+30;
  if(nphi<=-150) return 0;
  else if(nphi<=-90) return 1;
  else if(nphi<=-30) return 2;
  else if(nphi<=30) return 3;
  else if(nphi<=90) return 4;
  else if(nphi<=150) return 5;
  else return 0; 

}
