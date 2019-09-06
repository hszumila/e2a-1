#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring> 

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"

#include "event_Info.h"
#include "e2a_constants.h"

using namespace std;

const double deltaWidth = wD;
const TVector3 vBeam(0.,0.,4.461);
//Define some useful functions
double sq(double x){ return x*x; };

bool checkDelta(const double dMass);

double getMass(const TVector3 vN, const TVector3 vpi);

int main(int argc, char ** argv){

  if( argc != 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "leadCut /path/to/input/data/file /path/to/output/root/file \n";

    return -1;
    
  }
  
  //Get A trees. Open output file.
  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");

  //Make trees and histograms for data
  TTree * inTree = (TTree*)inputFile->Get("T");
  TTree * outTree = new TTree("T","Reduced Tree");

  cerr<<"Nucleus file has been opened from input file. \n";

  
  //Define varialbes needed for input tree and make tree
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

  cerr<<"input TTree opened from input file \n output TTree creadted. \n";

  //Define variable needed for output tree
  outTree->Branch("nParticles",&nPar,"nPar/I");
  outTree->Branch("Xb",&xB,"xB/D");
  outTree->Branch("Q2",&QSq,"QSq/D");           
  outTree->Branch("Part_type",parID,"parID[19]/I");
  outTree->Branch("mom_x",px,"px[19]/D");
  outTree->Branch("mom_y",py,"py[19]/D");
  outTree->Branch("mom_z",pz,"pz[19]/D");
  outTree->Branch("vtx_z_cor",vtxZCorr,"vtxZCorr[19]/D");
  
  TH1D * hist_Mass_binD[4];
  for(int i=0; i<4; i++){
      char temp[100];
      sprintf(temp,"hist_xB_Mass_D%d",i);
      hist_Mass_binD[i] = new TH1D(temp,"Delta;xB;Counts",40,1,1.4);
  }

  int fin = inTree->GetEntries();
  int numDeltas = 0;

  //Loop over TTree
  for(int i = 0; i < fin; i++){

    inTree->GetEntry(i);
    
    if(nPar > 19){
      cout<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }

    //Display completed
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }    

    //Now get the info so we can properly work with the array
    event_Info myInfo(nPar, parID, px, py, pz, vtxZCorr);
        
    //Loop over all particles and check for pions and nucleons
    int j = 1;
    int k = 1;
    while(j<myInfo.getNPar()){
      j++;
      while(k<myInfo.getNPar()){
	k++;

	//only do something if the particles give the correct combination
	if(!(myInfo.isNucleon(j) && myInfo.isPion(k))){
	  continue;
	}

	//Now get vector elements from the myInfo object
	TVector3 vN(myInfo.getPX(j),myInfo.getPY(j),myInfo.getPZ(j));
	TVector3 vpi(myInfo.getPX(k),myInfo.getPY(k),myInfo.getPZ(k));

	//Check if it is a delta, if it is, change the arrays  with
	//makeDelta. This will alter all of the arrays so now the
	//delta particle is at the back and all other particles have
	//been pushed forward in their place. Also after I find a delta,
	//I loop over the entire list again because now the ordering
	//has changed.
	double dMass = getMass(vN,vpi);
	if((checkDelta(dMass)) && myInfo.vtxMatch(j,k)){

	  if(myInfo.getDeltaType(j,k)==dmCode){hist_Mass_binD[0]->Fill(dMass);}
	  else if(myInfo.getDeltaType(j,k)==d0Code){hist_Mass_binD[1]->Fill(dMass);}
	  else if(myInfo.getDeltaType(j,k)==dpCode){hist_Mass_binD[2]->Fill(dMass);}
	  else if(myInfo.getDeltaType(j,k)==dppCode){hist_Mass_binD[3]->Fill(dMass);}

	  myInfo.addDelta(j,k);
	  j=0;
	  k=0;
	}
      }
    }
    

    //Now fill variables for writting out
    nPar=myInfo.getNPar();
    for(int l = 0; l < 19; l++){
      parID[l] = myInfo.getParID(l);
      px[l] = myInfo.getPX(l);
      py[l] = myInfo.getPY(l);
      pz[l] = myInfo.getPZ(l);
      vtxZCorr[l] = myInfo.getVTX(l);
    }

    numDeltas+=myInfo.getNumDeltas();
    if(myInfo.getNumDeltas() > 0){
    outTree->Fill();
    }
  }
    

  cout<<"I found: "<<numDeltas<<" deltas!\n";
  for(int l = 0; l < 4; l++){
    hist_Mass_binD[l]->Write();
  }
  outTree->Write();
  inputFile->Close();
  outputFile->Close();

  return 0;
  }


double getMass(const TVector3 vN, const TVector3 vpi){

  TVector3 vD_test = vN + vpi;
  double eN = sqrt(vN.Mag2()+sq(mN));
  double epi = sqrt(vpi.Mag2()+sq(mpc));
  double eD_test = eN + epi;
  double mD_test = sqrt(sq(eD_test)-vD_test.Mag2());
  
  return mD_test;    
}

bool checkDelta(const double dMass)
{
if((dMass > (mD-deltaWidth)) && (dMass < (mD+deltaWidth))){
      return true;
    }
  return false;

}

