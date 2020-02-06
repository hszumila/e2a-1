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
  cerr<< "Argumets: ./get_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Nucleus A for Ratio] [optional flags]\n\n"
      <<"Optional flags:\n"
      <<"-h: Help\n\n";
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
  int A = atoi(argv[3]);
  target_Info secondTargInfo(atoi(argv[4]));
  //Get target info
  target_Info targInfo(A);  
  bool verbose = true;
  
  int c;
  while ((c=getopt (argc-4, &argv[4], "h")) != -1) //First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_message();
	return -1;
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
  double finebinxB[] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.4,1.58,1.77,2};
  int finenumbinxB =  ( sizeof(finebinxB)/sizeof(finebinxB[0]) ) - 1;

  //Make a histogram list to make things easier
  vector<TH1*> hist_list;
  TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_xB =  new TH1D("hist_xB_Incl" ,"hist;xB;Counts",finenumbinxB,finebinxB);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",40,0,5);
  hist_list.push_back(hist_QSq);
  TH1D * hist_pMiss =  new TH1D("hist_pMiss" ,"hist;pMiss;Counts",40,0,2);
  hist_list.push_back(hist_QSq);

  TH1D * hist_eVTX =  new TH1D("hist_eVTX" ,"hist;eVTX;Counts",100,-10,10);
  hist_list.push_back(hist_eVTX);
  TH1D * hist_pVTX =  new TH1D("hist_pVTX" ,"hist;pVTX;Counts",100,-10,10);
  hist_list.push_back(hist_pVTX);


  TH1D * hist_xB_p[6];
  char temp[100];
  for(int i = 0; i<5; i++){
    sprintf(temp,"hist_xB_p_%d",i);
    hist_xB_p[i] = new TH1D(temp,"x_B Distribution;xB;counts",finenumbinxB,finebinxB);
    hist_list.push_back(hist_xB_p[i]);
  }
  sprintf(temp,"hist_xB_p_%d",5);
  hist_xB_p[5] = new TH1D(temp,"x_B Distribution;xB;counts",numbinxB,binxB);
  hist_list.push_back(hist_xB_p[5]);

  vector<TH2*> hist_list_2D;
  TH2D * mMiss_mX =  new TH2D("hist_mMiss_mX" ,"hist;mMiss;mX;Counts",100,0.8,1.6,400,0,4);
  hist_list_2D.push_back(mMiss_mX);  
  TH2D * xB_mMiss =  new TH2D("hist_xB_mMiss" ,"hist;xB;mMiss;Counts",finenumbinxB,finebinxB,100,0.8,1.6);
  hist_list_2D.push_back(xB_mMiss); 
  TH2D * xB_mX =  new TH2D("hist_xB_mX" ,"hist;xB;mX;Counts",finenumbinxB,finebinxB,400,0,4);
  hist_list_2D.push_back(xB_mX);  
  TH2D * hist_vtx = new TH2D("hist_VTX" ,"hist;eVTX;pVTX;Counts",100,-10,10,100,-10,10);
  hist_list_2D.push_back(hist_vtx);

  for(int i=0; i<hist_list_2D.size(); i++){
    hist_list_2D[i]->Sumw2();
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
    double weight = 1;
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);

    TVector3 ve = myInfo.getVector(0);
    TVector3 vq = vBeam - ve;
    double omega = vBeam.Mag() - ve.Mag();
    double phi = ve.Phi() * 180 / M_PI;
    double theta = ve.Theta() * 180 / M_PI;
    TVector3 vLead = myInfo.getVector(1);
    double thetaLead = vLead.Theta() * 180 / M_PI;
    TVector3 vMiss = vq - vLead;
    double massLead = mP;
    double eLead = sqrt( vLead.Mag2() + (massLead*massLead) );
    //double eMiss = massLead + omega - eLead;
    //double mX = (vMiss.Mag2() - (eMiss*eMiss))/(2 * eMiss);
    double mMiss = myInfo.getMassMiss(1);
    double eMiss = sqrt(vMiss.Mag2()+(mP*mP));
    double mX = eLead + eMiss - omega;
    double poq = myInfo.getPoQ(1);
    double thetapq = myInfo.getThetaPQ(1);
    double thetapMq = vq.Angle(vMiss) * (180/M_PI);
       
    
    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    

    
    //Only take particles that pass fiducials and acceptances for both targets
    if(myInfo.isProton(1)){
      if(targInfo.e_acc(ve) < accMin){continue;}
      if(targInfo.p_acc(vLead) < accMin){continue;}
      if(!targInfo.pass_semi_fid(ve,vLead)){ continue; }
      
      if(secondTargInfo.e_acc(ve) < accMin){continue;}
      if(secondTargInfo.p_acc(vLead) < accMin){continue;}
      if(!secondTargInfo.pass_semi_fid(ve,vLead)){ continue; }
    }



    weight = weight / targInfo.semi_acc(ve,vLead);
    weight = weight / targInfo.getTrans();
    weight = weight / (targInfo.getLum() * A);
    weight = weight / targInfo.getRadCorr(theta,xB);




    hist_eVTX->Fill(vtxZCorr[0],weight);
    hist_pVTX->Fill(vtxZCorr[1],weight);
    //if(!targInfo.vtxInRange(vtxZCorr[0],vtxZCorr[0],theta,thetaLead)){
    //  continue;
    //}
    hist_vtx->Fill(vtxZCorr[0],vtxZCorr[1],weight);

 
    
    hist_xB_p[0]->Fill(xB,weight);
    if(vMiss.Mag() < 0.1){continue;}
    hist_xB_p[1]->Fill(xB,weight);
    if(vMiss.Mag() < 0.13){continue;}
    hist_xB_p[2]->Fill(xB,weight);
    if(vMiss.Mag() < 0.16){continue;}
    hist_xB_p[3]->Fill(xB,weight);
    if(vMiss.Mag() < 0.2){continue;}
    hist_xB_p[4]->Fill(xB,weight);
    if(vMiss.Mag() < 0.3){continue;}
    hist_xB_p[5]->Fill(xB,weight);

    xB_mX->Fill(xB,mX,weight);
    mMiss_mX->Fill(mMiss,mX,weight);
    xB_mMiss->Fill(xB,mMiss,weight);
    hist_Cross->Fill(1,1);
    hist_xB->Fill(xB,weight);	  	  
    hist_QSq->Fill(QSq,weight);
    hist_pMiss->Fill(vMiss.Mag(),weight);
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
  for(int i=0; i<hist_list_2D.size(); i++){
    hist_list_2D[i]->Write();
  }
  outputFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}
