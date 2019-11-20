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

  TH1D * hist_mintpq =  new TH1D("hist_xB" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB);


  for(int k = 2; k < argc; k++){
    std::ifstream inFile(argv[k]);
    string temp;
    double mintpq, maxtpq, minpoq, minpmiss, maxpmiss, minmass, maxmass;
    inFile >> temp >> temp >> temp >> temp >> temp >> temp;
    inFile >> maxtpq >> minpoq >> minpmiss >> maxpmiss >> minmass >> maxmass;

  }
  outputFile->Close();
  return 0;
}
