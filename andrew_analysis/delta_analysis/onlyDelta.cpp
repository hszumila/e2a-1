#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring> 
#include <unistd.h>

#include "TMath.h"
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
  cerr<< "Argumets: ./leadCut /path/to/input/skim/file /path/to/output/root/file\n\n"
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
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    }
  
  //Read in arguments.
  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");


  bool verbose = false;
   
  int c;
  while ((c=getopt (argc-2, &argv[2], "hv")) != -1) //First two arguments are not optional flags.
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


  //Make trees and histograms for data
  TTree * inTree = (TTree*)inputFile->Get("T");
  TTree * outTree = new TTree("T","Reduced Tree");

  vector<TH1*> hist_list;
  TH1D * hist_Mass_m =  new TH1D("hist_Mass_m" ,"hist;Mass_m;Counts",100,1,1.9);
  hist_list.push_back(hist_Mass_m);
  TH1D * hist_Mass_0 =  new TH1D("hist_Mass_0" ,"hist;Mass_0;Counts",100,1,1.9);
  hist_list.push_back(hist_Mass_0);
  TH1D * hist_Mass_p =  new TH1D("hist_Mass_p" ,"hist;Mass_p;Counts",100,1,1.9);
  hist_list.push_back(hist_Mass_p);
  TH1D * hist_Mass_pp =  new TH1D("hist_Mass_pp" ,"hist;Mass_pp;Counts",100,1,1.9);
  hist_list.push_back(hist_Mass_pp);


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
 int numDelm = 0;
 int numDel0 = 0;
 int numDelp = 0;
 int numDelpp = 0;
  
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
    double delMass = myInfo.findAndMergeDeltas();
    if(delMass < 0){ continue; }
    myInfo.clearNonDelta();

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
    //Now apply some counting
    if(myInfo.getParID(1) == dmCode){
      numDelm++;
      hist_Mass_m->Fill(delMass);
    }
    if(myInfo.getParID(1) == d0Code){
      numDel0++;
      hist_Mass_0->Fill(delMass);
    }
    if(myInfo.getParID(1) == dpCode){
      numDelp++;      
      hist_Mass_p->Fill(delMass);
    }
    if(myInfo.getParID(1) == dppCode){
      numDelpp++;
      hist_Mass_pp->Fill(delMass);
    }
    outTree->Fill();
    nEvents++;
  }

  cout<<"Finished filling tree for output"<<endl;
  cout<<"Number of Events with Deltas: "<<nEvents<<endl;
  cout<<"Number of Events with Delta Minus: "<<numDelm<<endl;
  cout<<"Number of Events with Delta Naught: "<<numDel0<<endl;
  cout<<"Number of Events with Delta Plus: "<<numDelp<<endl;
  cout<<"Number of Events with Delta Plus Plus: "<<numDelpp<<endl;
  outTree->Write();
  for(int k=0; k<hist_list.size(); k++){
    hist_list[k]->Write();
  }
  inputFile->Close();
  outputFile->Close();

  return 0;
}
  
    //if(mergeD){myInfo.findAndMergeDeltas();}
    //if(mergeD && (myInfo.getNumDeltas() < 1)){continue;}
    
    //int leadIndex = myInfo.getWhichLead();
    
    //if(checkLeadD && !myInfo.doesNucleonMatchPion(leadIndex)){continue;}
    
    //if(leadIndex > 0){
    //  myInfo.setLead(leadIndex);
