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
int getPMissBin(double x){
  if(x<0.05){return 0;}
  if(x<0.1){return 1;}
  if(x<0.2){return 2;}
  if(x<0.4){return 3;}
  if(x<0.6){return 4;}
  else return -1;
}
int getaMissBin(double x){
  if(x<0.7){return 0;}
  if(x<1){return 1;}
  if(x<1.3){return 2;}
  else return 3;
}

int getPTCut(double x){
  if(x<0.05){return 0;}
  if(x<0.1){return 1;}
  if(x<0.3){return 2;}
  else return 3;
}
int getQCut(double x){
  if(x<1.0){return -1;}
  if(x<1.1){return 0;}
  if(x<1.2){return 1;}
  if(x<1.3){return 2;}
  if(x<1.4){return 3;}
  if(x<1.5){return 4;}
  else return 5;
}
int getPMissCut(double x){ 
  if(x<0.05){return 0;}
  if(x<0.1){return 1;}
  if(x<0.15){return 2;}
  if(x<0.2){return 3;}
  if(x<0.25){return 4;}
  if(x<0.3){return 5;}
  else return 6;
}

void help_message()
{
  cerr<< "Argumets: ./get_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Nucleus A for Ratio] [optional flags]\n\n"
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
  double binxB[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.58,1.77,2};
  //double binxB[] = {0,0.3,0.5,0.7,0.9,1.1,1.3,1.58,2};
  int numbinxB = ( sizeof(binxB)/sizeof(binxB[0]) ) - 1;
  double finebinxB[] = {0.4,0.5,0.6,0.68,0.74,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.4,1.58,1.77,2};
  int finenumbinxB =  ( sizeof(finebinxB)/sizeof(finebinxB[0]) ) - 1;
  double binalpha[] = {0.4,0.5,0.6,0.68,0.74,0.8,0.85,0.88,0.91,0.94,0.97,1,1.03,1.06,1.09,1.12,1.15,1.2,1.26,1.32,1.4,1.5,1.6};
  int numbinalpha = ( sizeof(binalpha)/sizeof(binalpha[0]) ) - 1;

  //Make a histogram list to make things easier
  vector<TH1*> hist_list;
  TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_xB =  new TH1D("hist_xB" ,"hist;xB;Counts",finenumbinxB,finebinxB);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",40,0,5);
  hist_list.push_back(hist_QSq);
  TH1D * hist_pLead =  new TH1D("hist_pLead" ,"hist;pLead;Counts",40,0,5);
  hist_list.push_back(hist_pLead);
  TH1D * hist_tLead =  new TH1D("hist_tLead" ,"hist;tLead;Counts",90,0,180);
  hist_list.push_back(hist_tLead);
  TH1D * hist_pMiss =  new TH1D("hist_pMiss" ,"hist;pMiss;Counts",24,0,0.6);
  hist_list.push_back(hist_pMiss);

  TH1D * hist_eVTX =  new TH1D("hist_eVTX" ,"hist;eVTX;Counts",100,-10,10);
  hist_list.push_back(hist_eVTX);
  TH1D * hist_pVTX =  new TH1D("hist_pVTX" ,"hist;pVTX;Counts",100,-10,10);
  hist_list.push_back(hist_pVTX);


  TH1D * hist_xB_pCut[7];
  TH1D * hist_xB_pBin[5];
  char temp[100];
  for(int i = 0; i<5; i++){
    sprintf(temp,"hist_xB_pCut_%d",i);
    hist_xB_pCut[i] = new TH1D(temp,"x_B Distribution;xB;counts",finenumbinxB,finebinxB);
    hist_list.push_back(hist_xB_pCut[i]);

    sprintf(temp,"hist_xB_pBin_%d",i);
    hist_xB_pBin[i] = new TH1D(temp,"x_B Distribution;xB;counts",finenumbinxB,finebinxB);
    hist_list.push_back(hist_xB_pBin[i]);
  }

  for(int i = 5; i<7; i++){
    sprintf(temp,"hist_xB_pCut_%d",i);
    hist_xB_pCut[i] = new TH1D(temp,"x_B Distribution;xB;counts",numbinxB,binxB);
    hist_list.push_back(hist_xB_pCut[i]);
  }

  TH1D * hist_pp_pCut[7];
  TH1D * hist_pt_pCut[7];
  TH1D * hist_aMiss[7];  
  for(int i = 0; i<7; i++){
    sprintf(temp,"hist_pMiss_parallel_pCut_%d",i);
    hist_pp_pCut[i] = new TH1D(temp,"pMiss parallel Distribution;pp;counts",40,-0.6,0.6);
    hist_list.push_back(hist_pp_pCut[i]);

    sprintf(temp,"hist_pMiss_transverse_pCut_%d",i);
    hist_pt_pCut[i] = new TH1D(temp,"pMiss transverse Distribution;pt;counts",20,0,0.6);
    hist_list.push_back(hist_pt_pCut[i]);

    sprintf(temp,"hist_aMiss_pCut_%d",i);
    hist_aMiss[i] =  new TH1D(temp ,"hist;aMiss;Counts",numbinalpha,binalpha);
    hist_list.push_back(hist_aMiss[i]);

  }

  TH1D * hist_pt_aBin[4];
  for(int i = 0; i<7; i++){
    sprintf(temp,"hist_pt_aMissBin_%d",i);
    hist_pt_aBin[i] = new TH1D(temp,"pLead transverse Distribution;pt;counts",20,0,0.6);
    hist_list.push_back(hist_pt_aBin[i]);
  }


  TH1D * hist_xB_qCut[6];
  for(int i = 0; i<6; i++){
    sprintf(temp,"hist_xB_qCut_%d",i);
    hist_xB_qCut[i] = new TH1D(temp,"x_B Distribution;xB;counts",numbinxB,binxB);
    hist_list.push_back(hist_xB_qCut[i]);
  }

  TH1D * hist_pMiss_ptCut[4];
  TH1D * hist_xB_ptCut[4];
  TH1D * hist_aMiss_ptCut[4];
  for(int i = 0; i<4; i++){
    sprintf(temp,"hist_pMiss_ptCut_%d",i);
    hist_pMiss_ptCut[i] = new TH1D(temp,"p_Miss Distribution;p_Miss;counts",24,0,0.6);
    hist_list.push_back(hist_pMiss_ptCut[i]);

    sprintf(temp,"hist_xB_ptCut_%d",i);
    hist_xB_ptCut[i] = new TH1D(temp,"x_B Distribution;x_B;counts",numbinxB,binxB);
    hist_list.push_back(hist_xB_ptCut[i]);

    sprintf(temp,"hist_aMiss_ptCut_%d",i);
    hist_aMiss_ptCut[i] = new TH1D(temp,"x_B Distribution;x_B;counts",numbinalpha,binalpha);
    hist_list.push_back(hist_aMiss_ptCut[i]);
  }




  vector<TH2*> hist_list_2D;
  TH2D * xB_pMiss =  new TH2D("hist_xB_pMiss" ,"hist;xB;pMiss;Counts",finenumbinxB,finebinxB,100,0,2);
  hist_list_2D.push_back(xB_pMiss); 
  TH2D * xB_mMiss =  new TH2D("hist_xB_mMiss" ,"hist;xB;mMiss;Counts",finenumbinxB,finebinxB,100,0.8,1.6);
  hist_list_2D.push_back(xB_mMiss); 
  TH2D * xB_QSq =  new TH2D("hist_xB_QSq" ,"hist;xB;QSq;Counts",finenumbinxB,finebinxB,100,0,5);
  hist_list_2D.push_back(xB_QSq); 
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
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);

    TVector3 ve = myInfo.getVector(0);
    TVector3 vq = vBeam - ve;
    double omega = vBeam.Mag() - ve.Mag();
    double phi = ve.Phi() * 180 / M_PI;
    double eventSector = targInfo.getSecPhi(phi);
    double theta = ve.Theta() * 180 / M_PI;
    TVector3 vLead = myInfo.getVector(1);
    double pLead = vLead.Mag();
    double thetaLead = vLead.Theta() * 180 / M_PI;
    TVector3 vMiss = vq - vLead;
    double massLead = mP;
    double eLead = sqrt( vLead.Mag2() + (massLead*massLead) );
    double mMiss = myInfo.getMassMiss(1);
    double poq = myInfo.getPoQ(1);
    double thetapq = myInfo.getThetaPQ(1);
    double thetapMq = vq.Angle(vMiss) * (180/M_PI);
    double pMiss_p = vMiss.Dot(vq.Unit());
    TVector3 pLead_t_vec = vLead.Perp(vq);
    TVector3 pMiss_t_vec = vMiss.Perp(vq);
    double pMiss_t = pMiss_t_vec.Mag();
    double pLead_t = pLead_t_vec.Mag();

    double mA = targInfo.getMass();
    TVector3 qUnit = vq.Unit();
    double aLead = (eLead - vLead.Dot(qUnit)) / (mA/(double)A);
    double aq = (omega - vq.Dot(qUnit)) / (mA/(double)A);
    double aMiss = aLead - aq;
    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cout << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    

    
    //Only take particles that pass fiducials and acceptances for both targets
    //    if(myInfo.isProton(1)){
    if(targInfo.e_acc(ve) < accMin){continue;}
    if(targInfo.p_acc(vLead) < accMin){continue;}
    if(!targInfo.pass_semi_fid(ve,vLead)){ continue; }
      
    if(secondTargInfo.e_acc(ve) < accMin){continue;}
    if(secondTargInfo.p_acc(vLead) < accMin){continue;}
    if(!secondTargInfo.pass_semi_fid(ve,vLead)){ continue; }
    //    }
    //if(doOneSec && (eventSector != mySec)){ continue; }
      
    //Establish vertex Cuts
    hist_eVTX->Fill(vtxZCorr[0],1);
    hist_pVTX->Fill(vtxZCorr[1],1);
    hist_vtx->Fill(vtxZCorr[0],vtxZCorr[1],1);
    //if(!targInfo.eVTXInRange(vtxZCorr[0])){
    if(!targInfo.semiFixedVTXInRange(myInfo,1)){
      continue;
    }


    //Apply the corrections to the cross sections
    // Acceptance correction
    // Transparancy correcton
    // Radiative correction
    // Normalize with luminosity and A
    double weight = 1;
    weight = weight / targInfo.semi_acc(ve,vLead);
    weight = weight / targInfo.getTrans();
    weight = weight / targInfo.getRadCorr(theta,xB);
    weight = weight / (targInfo.getLum() * A);

    counter++;
    //They all get these cuts

    /*    if(!((weight>0) && (weight <10000000))){
      cout<<"xB = " << xB << "\n";
      cout<<"weight = " << weight << "\n";
      cout<<"acceptance = " << targInfo.semi_acc(ve,vLead) << "\n";
      cout<<"transparency = " << targInfo.getTrans() << "\n";
      cout<<"rad = " << targInfo.getRadCorr(theta,xB) << "\n";
      cout<<"lum = " << (targInfo.getLum() * A) << "\n";
      }*/

    if(vMiss.Mag() > 0.6){continue;}

    if(vMiss.Mag()>0.3){
      int QStart = getQCut(QSq);
      if(QStart >-1){
	for(int l = 0; l <= getQCut(QSq); l++){
	hist_xB_qCut[l]->Fill(xB,weight);
	}
      }
    }

    if(QSq<1.4){continue; }
    //No pMiss cut
    hist_xB->Fill(xB,weight);	  	  
    xB_pMiss->Fill(xB,vMiss.Mag(),weight);
    hist_pMiss->Fill(vMiss.Mag(),weight);

    //Bins in pMiss and aMiss
    hist_pt_aBin[getaMissBin(aMiss)]->Fill(pLead_t);
    hist_xB_pBin[getPMissBin(vMiss.Mag())]->Fill(xB,weight);

    //Do a pt cut in lieu of a pMiss cut
    for(int m = getPTCut(pLead_t); m <= 3; m++){
      hist_pMiss_ptCut[m]->Fill(vMiss.Mag(),weight);
      hist_xB_ptCut[m]->Fill(xB,weight);
      hist_aMiss_ptCut[m]->Fill(aMiss,weight);
    }

    //Contant pMiss cut
    if(vMiss.Mag()>0.3){

      hist_QSq->Fill(QSq,weight);
      hist_pLead->Fill(pLead,weight);
      hist_tLead->Fill(thetaLead,weight);

      xB_mMiss->Fill(xB,mMiss,weight);
      xB_QSq->Fill(xB,QSq,weight);
      hist_Cross->Fill(1,weight);

    }

    //Apply variable pMiss Cuts
    for(int l = 0; l <= getPMissCut(vMiss.Mag()); l++){
      hist_xB_pCut[l]->Fill(xB,weight);
      hist_pp_pCut[l]->Fill(pMiss_p,weight);
      hist_pt_pCut[l]->Fill(pMiss_t,weight);
      hist_aMiss[l]->Fill(aMiss,weight);
    }

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
