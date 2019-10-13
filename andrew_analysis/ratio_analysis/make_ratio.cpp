#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"

#include "target_Info.h"

using namespace std;

double sq(double x){
  return x*x;
}

int main(int argc, char ** argv){

  if( argc != 4){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "make_ratio /path/to/input/numerator/file /path/to/input/denominator/file /path/to/output/file \n";

    return -1;
    
  }
  
  //Get D and A trees. Open output file.
  TFile * inputNumFile = new TFile(argv[1]);
  TList * numList = (TList*)inputNumFile->GetListOfKeys();
  TFile * inputDenFile = new TFile(argv[2]);
  TList * denList = (TList*)inputDenFile->GetListOfKeys();
  TFile * outputFile = new TFile(argv[3],"RECREATE");
  cerr<<"Files have been opened\n";
  
  outputFile->cd();
  //Loop over objects in the numerator
  for(int i = 0; i < numList->GetEntries(); i++){
    //Loop over objects in the denominator
    for(int j = 0; j < denList->GetEntries(); j++){


      //Now get the names of the objects, check what types of objects they are, and check if they are compatible
      const char * sN=numList->At(i)->GetName();
      const char * sD=denList->At(j)->GetName();
      bool areTH1 =( (inputNumFile->Get(sN)->InheritsFrom(TH1D::Class())) && (inputDenFile->Get(sD)->InheritsFrom(TH1D::Class())) ); 
      bool areTH2 =( (inputNumFile->Get(sN)->InheritsFrom(TH2D::Class())) && (inputDenFile->Get(sD)->InheritsFrom(TH2D::Class())) ); 
      if(!(strcmp(sN,sD)==0)){
	continue;
      }
    
      //If they are nD histos, then divide using nD histos
      if(areTH1){
	cout<<sN<<"\n";
	//Get Histograms
	TH1D * num_hist = (TH1D*)inputNumFile->Get(sN);
	TH1D * den_hist = (TH1D*)inputDenFile->Get(sD);
	TH1D * ratio_hist = (TH1D*)num_hist->Clone(sD);
	ratio_hist->Divide(den_hist);
	//	cout<<"Num Error: "<<sq((num_hist->GetBinError(6))/(num_hist->GetBinContent(6)))<<"     Den Error: "<<sq((den_hist->GetBinError(6))/(den_hist->GetBinContent(6)))<<"     Ratio Error: "<<sq((ratio_hist->GetBinError(6))/(ratio_hist->GetBinContent(6)))<<"\n";
	ratio_hist->Write();
      }    
      else if (areTH2){
	cout<<sN<<"\n";
	//Get Histograms
	TH2D * num_hist = (TH2D*)inputNumFile->Get(sN);
	TH2D * den_hist = (TH2D*)inputDenFile->Get(sD);
	TH2D * ratio_hist = (TH2D*)num_hist->Clone(sD);
	ratio_hist->Divide(den_hist);
	ratio_hist->Write();	
      }

    }
  }
    
  cerr<<"Ratio made\n";
 
  //Now use ttrees to get the total cross section
  TH1D * numHist = (TH1D*)inputNumFile->Get("totalCross");
  double numCross = numHist->GetBinContent(1);
  TH1D * denHist = (TH1D*)inputDenFile->Get("totalCross");
  double denCross = denHist->GetBinContent(1);
  cout<<"The cross section for this ratio is: "<<(numCross/denCross)<<"\n";
  
  //Close Files
  outputFile->Close();
  inputNumFile->cd();
  inputNumFile->Close();
  inputDenFile->cd();
  inputDenFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
  /* //Code for using a TGraphAsymmErrors
  num_xB->Scale(0.01);  
  TGraphAsymmErrors * ratio_xB = new TGraphAsymmErrors();
  //ratio_xB->Divide(num_xB,den_xB);
  ratio_xB->Divide(num_xB,den_xB,"cl=0.683 b(1,1) mode");
  */  
