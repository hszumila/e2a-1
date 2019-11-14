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
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "target_Info.h"

using namespace std;

double sq(double x){
  return x*x;
}

int main(int argc, char ** argv){

  if( argc < 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "make_ratio /path/to/output/file /path/to/all/input/files \n";

    return -1;
    
  }

  TFile * outputFile = new TFile(argv[1],"RECREATE");
  //Counter for random names
  int ctr = 0;
  TFile * File1 = new TFile(argv[2]);
  TList * List1 = (TList*)File1->GetListOfKeys();

  //Loop over all of the objects in the first file to make a tgraph of those objects
  for(int i = 0; i < List1->GetEntries(); i++){
    const char * Name1 = List1->At(i)->GetName();
    //Only use it if it is a TH1D
    if(!(File1->Get(Name1)->InheritsFrom(TH1D::Class()))){continue;}
    TH1D * Hist1 = (TH1D*)File1->Get(Name1);
    cout<<"Working on: "<<Name1<<"\n";
    
    //If it is a TH1D then make a set of histograms for it so you can find the mean and stuff
    int numBin = Hist1->GetNbinsX();
    TH1D * histOfPoints[numBin];
    for(int j = 0; j < numBin; j++){
      char temp[100];
      sprintf(temp,"Hist_%i",ctr);
      ctr++;
      histOfPoints[j] = new TH1D(temp,temp,120,0,4);
    }
    
    //Now fill the histogram with points
    //First loop through each file
    for(int k = 2; k < argc; k++){
      TFile * inFile = new TFile(argv[k]);
      TList * inList = (TList*)inFile->GetListOfKeys();
      cout<<"Looking at File: " << argv[k] << "\n";
	
      //Then loop through each object
      for(int l = 0; l < inList->GetEntries(); l++){
	const char * inName=inList->At(l)->GetName();
	//Only use it if it is a TH1D
	if(inFile->Get(inName)->InheritsFrom(TH1D::Class())){
	  //Only use it if it matches the TH1D we are trying to get
	  if(strcmp(Name1,inName)==0){
	    TH1D * inHist = (TH1D*)inFile->Get(inName);
	    //Only use if the number of bins matches
	    if((inHist->GetNbinsX()) == numBin){
	      //Now loop over each point
	      for(int m = 0; m < numBin; m++){
		histOfPoints[m]->Fill(inHist->GetBinContent(m+1));
	      }
	
	    }

	  }

	}
      }

      inFile->Close();
    }

    TGraphAsymmErrors * histGraph = new TGraphAsymmErrors();
    histGraph->SetName(Name1);
    //Then loop through each bin to get the mean and std
    for(int n = 0; n < numBin; n++){
      double xValue, mean, sigma;
      xValue = Hist1->GetBinCenter(n+1);
      mean = histOfPoints[n]->GetMean();
      sigma = histOfPoints[n]->GetStdDev();
      //Make the fit
      TF1 * myfit = new TF1(Form("Fit_%i",ctr),"gaus",0,4);
      ctr++;
      myfit->SetParameter(0,histOfPoints[n]->GetMaximum());
      myfit->SetParameter(1,mean);
      myfit->SetParameter(2,sigma);
      TFitResultPtr fitPoint = histOfPoints[n]->Fit(myfit,"qeSrn","",0,4);
      int fitStatus = fitPoint;
      //Change only if fit passes
      if(fitStatus == 0){
	mean = fitPoint->Parameter(1);
	sigma = fitPoint->Parameter(2);
      }
      else{
	cout<<"Failed Fit...\n";
      }
      //Fill tgraphs with the results of the fit
      histGraph->SetPoint(histGraph->GetN(),xValue,mean);
      histGraph->SetPointError(histGraph->GetN()-1,0,0,sigma,sigma);
    }

    //Now write out the tgraphs
    outputFile->cd();
    for(int p = 0; p < numBin; p++){
      histOfPoints[p]->Write();
    }
    histGraph->Write();
    
  }

  outputFile->Close();
  return 0;
}
