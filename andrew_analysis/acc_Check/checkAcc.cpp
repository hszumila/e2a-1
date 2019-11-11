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
const double minPe = 1.07;
const double maxPe = 4.461;
const double minTheta = 14;
const double maxTheta = 55;
const int binSize = 500;
double sq(double x){
  return x*x;
}
int getSec(const double phi);

void help_message()
{
  cerr<< "Argumets: ./incl_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-v: Verbose\n"
      <<"-n: Set a custom minimum to the Z vertex [cm]\n"
      <<"-x: Set a custom maximum to the Z vertex [cm]\n"
      <<"-m: Do not apply maps to weight\n\n";
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
  TFile * outputFile = new TFile(argv[2],"RECREATE");
  int A = atoi(argv[3]);
  //Get target info
  target_Info targInfo(A);
  
  bool verbose = false;
  bool doMaps = true;
  
  int c;
  while ((c=getopt (argc-3, &argv[3], "hvn:x:m")) != -1) //First two arguments are not optional flags.
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
      case 'm':
	doMaps = false;
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cerr<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");

  //Make a histogram list to make things easier
  TH1D * hist_eAcc =  new TH1D("hist_eAcc" ,"hist;eAcc;Counts",100,0,1);
  hist_eAcc->Sumw2();

  vector<TH2*> hist_list;
  TH2D * pvt_exp[6];
  TH2D * pvt_bad[6];  
  TH2D * pvt_map[6];
  TH2D * tvf_exp[6];
  TH2D * tvf_bad[6];  
  TH2D * tvf_map[6];


  for(int i=0; i<6; i++){
      char temp[100];

      sprintf(temp,"pe_v_theta_exp_%d",i);
      pvt_exp[i] = new TH2D(temp,"pe_v_theta;pe;theta;Counts",100,minPe,maxPe,100,minTheta,maxTheta);
      hist_list.push_back(pvt_exp[i]);

      sprintf(temp,"pe_v_theta_bad_%d",i);
      pvt_bad[i] = new TH2D(temp,"pe_v_theta;pe;theta;Counts",100,minPe,maxPe,100,minTheta,maxTheta);
      hist_list.push_back(pvt_bad[i]);

      sprintf(temp,"pe_v_theta_map_%d",i);
      pvt_map[i] = new TH2D(temp,"pe_v_theta;pe;theta;Counts",100,minPe,maxPe,100,minTheta,maxTheta);
      hist_list.push_back(pvt_map[i]);

      sprintf(temp,"theta_v_phi_exp_%d",i);
      tvf_exp[i] = new TH2D(temp,"theta_v_phi;pe;theta;Counts",100,minTheta,maxTheta,1000,-180,180);
      hist_list.push_back(tvf_exp[i]);

      sprintf(temp,"theta_v_phi_bad_%d",i);
      tvf_bad[i] = new TH2D(temp,"theta_v_phi;pe;theta;Counts",100,minTheta,maxTheta,1000,-180,180);
      hist_list.push_back(tvf_bad[i]);

      sprintf(temp,"theta_v_phi_map_%d",i);
      tvf_map[i] = new TH2D(temp,"theta_v_phi;pe;theta;Counts",100,minTheta,maxTheta,100,-180,180);
      hist_list.push_back(tvf_map[i]);
      
  }   
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

    pvt_exp[sec]->Fill(pe,theta,1);
    tvf_exp[sec]->Fill(theta,phi,1);
    hist_eAcc->Fill(targInfo.e_acc(ve));

    if(targInfo.e_acc(ve)<0.1){
      pvt_bad[sec]->Fill(pe,theta,1);
      tvf_bad[sec]->Fill(theta,phi,1);
      //      cout<<pe<<" "<<theta<<" "<<phi<<" \n";
    }
  }
  

  //Now scan over all space to get the map
  for(int i = 0; i < binSize; i++){
    //Display completed
    if(verbose && (((i*100)%binSize)==0)){
      cerr << (i*100/binSize) <<"% complete \n";
    }

    for(int j = 0; j < binSize; j++){
      for(int k = 0; k < 100; k++){
	double pe = minPe + (((maxPe-minPe)/binSize) * ((double)i + 0.5));
	double theta = minTheta + (((maxTheta-minTheta)/binSize) * ((double)j) + 0.5);
	double phi = -180 + ((360/100)*((double)k + 0.5));
	int sec = getSec(phi);
	
	TVector3 ve;//(pe*sin(theta)*cos(phi), pe*sin(theta)*sin(phi), pe*cos(theta));
	ve.SetMagThetaPhi(pe,(theta * M_PI / 180),(phi * M_PI / 180));
	double acc=targInfo.e_acc(ve);	
	pvt_map[sec]->Fill(pe,theta,acc);
	tvf_map[sec]->Fill(theta,phi,acc);

      }
    }
  }
  cerr<<"Finished filling histogram\n";

   
  inputFile->cd();
  inputFile->Close();
  //Now write out
  outputFile->cd();
  hist_eAcc->Write();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  outputFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
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
