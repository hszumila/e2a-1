#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring> 
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"

#include "event_Info.h"
#include "srcCut_Info.h"
#include "e2a_constants.h"

using namespace std;

double sq(double x){ return x*x; };

void help_message()
{
  cerr<< "Argumets: ./apply_SRC_Cut /path/to/input/skim/file /path/to/output/root/file\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-v: Verbose\n"
      <<"-x: Set minimum xB\n"
      <<"-p: Set minimum pMiss [GeV]\n\n";
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
  
  //Read in arguments.
  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");
  srcCut_Info myCut;
  myCut.setOnlyLeadProton();
  myCut.makeEG2Cut();  

  bool verbose = false;
   
  int c;
  while ((c=getopt (argc-4, &argv[4], "hvp:x:")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
      case 'v':
	verbose = true;
	break;
      case 'p':
	myCut.setMinPMissCut(atof(optarg));
	break;
      case 'x':
	myCut.setMinXBCut(atof(optarg));
	break;
      case '?':
	return -1;
      default:
	abort();
      }


  //Make trees and histograms for data
  TTree * inTree = (TTree*)inputFile->Get("T");
  TTree * outTree = new TTree("T","Reduced Tree");

  cerr<<"Nucleus file has been opened from input file. \n";

  
  //Define varialbes needed for input tree and make tree
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

  //Define branches for output tree
  outTree->Branch("nParticles",&nPar,"nPar/I");  
  outTree->Branch("Xb",&xB,"xB/D");            
  outTree->Branch("Q2",&QSq,"QSq/D");           
  outTree->Branch("Part_type",parID,"parID[19]/I");   
  outTree->Branch("mom_x",momx,"momx[19]/D");          
  outTree->Branch("mom_y",momy,"momy[19]/D");          
  outTree->Branch("mom_z",momz,"momz[19]/D");          
  outTree->Branch("vtx_z_cor",vtxZCorr,"vtxZCorr[19]/D");
  

  cerr<<"input TTree opened from input file \n output TTree creadted. \n";
  
  int fin = inTree->GetEntries();
  int nEvents = 0;
  int numLeadP = 0;
  int numLeadN = 0;
  int numDouble = 0;
  
  //Loop over TTree
  for(int i = 0; i < fin; i++){

    inTree->GetEntry(i);
    
    //Display completed
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }    
    
    if(nPar > 19){
      cout<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }

    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);
    int passCut = myCut.passCutReorder(myInfo);
    if(passCut == 0){ continue; }
    if(passCut == 2){
      numDouble++;
      continue;
    }

    //Now fill variables that have changed for writting out
    nPar=myInfo.getNPar();
    xB=myInfo.getXB();
    QSq=myInfo.getQSq();
    for(int l = 0; l < nPar; l++){
      parID[l] = myInfo.getParID(l);
      momx[l] = myInfo.getPX(l);
      momy[l] = myInfo.getPY(l);
      momz[l] = myInfo.getPZ(l);
      vtxZCorr[l] = myInfo.getVTX(l);
    }
    //Now now count things up to make sure it all works
    if(myInfo.isProton(1)){numLeadP++;}
    if(myInfo.isNeutron(1)){numLeadN++;}
    outTree->Fill();
    nEvents++;
  }

  cout<<"Finished filling tree for output"<<endl;
  cout<<"Number of Events with lead Nucleons: "<<nEvents<<endl;
  cout<<"Number of Events with Lead Proton: "<<numLeadP<<endl;
  cout<<"Number of Events with Lead Neutron: "<<numLeadN<<endl;
  cout<<"Number of Events with more than one Lead Nucleon: "<<numDouble<<endl;
  outTree->Write();
  inputFile->Close();
  outputFile->Close();

  return 0;
}
