#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 
#include <unistd.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"

#include "target_Info.h"
#include "event_Info.h"

using namespace std;

//Define some globale functions and constants
double accMin = 0.4;
double sq(double x){
  return x*x;
}
int getPMin(const double PMiss);
int getPMax(const double PMiss);
int getQSqBin(const double QSq);
int getXBBin(const double XB);
int getPMissBin(const double PMiss);
int getSec(const double phi);

void help_message()
{
  cerr<< "Argumets: ./get_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Nucleus A for Ratio] [Data Type] [optional flags]\n\n"
      <<"Data Types:\n"
      <<"1: Inclusive (only electrons)\n"
      <<"2: Semi-Inclusive (electron and lead nucleon)\n"
      <<"4: Semi-Inclusive with Delta\n\n"      
      <<"Optional flags:\n"
      <<"-h: Help\n"
      <<"-v: Verbose\n"
      <<"-n: Set a custom minimum to the Z vertex [cm]\n"
      <<"-x: Set a custom maximum to the Z vertex [cm]\n"
      <<"-s: Look at one sector [0-5]\n"
      <<"-Q: Set minimum QSq value [GeV] (default = 1.4)\n"
      <<"-t: Do not apply transparency factors for the lead proton\n"
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
  if (argc < 6)
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
  int TypeData = atoi(argv[5]);
  //Get target info
  target_Info targInfo(A);  
  bool verbose = false;
  bool doMaps = true;
  bool doTrans = true;
  bool oneSec = false;
  double minQSq = 1.4;
  int secChoice = -1;
  
  int c;
  while ((c=getopt (argc-5, &argv[5], "hvn:x:ms:tQ:")) != -1) //First two arguments are not optional flags.
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
      case 's':
	oneSec = true;
	secChoice = atof(optarg);
	break;
      case 't':
	doTrans = false;
	break;
      case 'Q':
	minQSq = atof(optarg);
	break;
      case '?':
	return -1;
      default:
	abort();
      }
  
  
  cerr<<"Files have been opened\n";

  //Make Trees and histograms
  TTree * inTree = (TTree*)inputFile->Get("T");
  double binxBIncl[] = {0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.38,1.58,1.79,2};
  int numbinxBIncl = 22;
  double binxB[] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.38,1.58,1.79,2};
  int numbinxB = 30;

  //Make a histogram list to make things easier
  vector<TH1*> hist_list;
  vector<TH2*> twoD_list;
  TH2D * poq_thetapq_0 = new TH2D("poq_thetapq_0","poq_thetapq_0;poq;thetapq;Counts",24,0,1.2,60,0,60);
  twoD_list.push_back(poq_thetapq_0);
  TH2D * poq_thetapq_1 = new TH2D("poq_thetapq_1","poq_thetapq_1;poq;thetapq;Counts",24,0,1.2,60,0,60);
  twoD_list.push_back(poq_thetapq_1);
  TH2D * xB_mMiss_All[6][3];
  TH2D * xB_mMiss_cutPoQ[6][3];
  TH2D * xB_poq_All[6][3];
  TH2D * xB_poq_cutMMiss[6][3];
  TH1D * hist_xB_p[6][3];
  for(int i=0; i<6; i++){
    for(int j=0; j<3; j++){
      char temp[100];

      sprintf(temp,"xB_mMiss_All_%d_%d",i,j);
      xB_mMiss_All[i][j] = new TH2D(temp,"x_B vs m_Miss;xB;mMiss",numbinxB,binxB,40,0.5,2.5);
      twoD_list.push_back(xB_mMiss_All[i][j]);

      sprintf(temp,"xB_mMiss_cutPoQ_%d_%d",i,j);
      xB_mMiss_cutPoQ[i][j] = new TH2D(temp,"x_B vs m_Miss with PoQ Cut;xB;mMiss",numbinxB,binxB,40,0.5,2.5);
      twoD_list.push_back(xB_mMiss_cutPoQ[i][j]);

      sprintf(temp,"xB_poq_All_%d_%d",i,j);
      xB_poq_All[i][j] = new TH2D(temp,"x_B vs p/q;xB;p/q",numbinxB,binxB,20,0.4,1.4);
      twoD_list.push_back(xB_poq_All[i][j]);

      sprintf(temp,"xB_poq_cutMMiss_%d_%d",i,j);
      xB_poq_cutMMiss[i][j] = new TH2D(temp,"x_B vs p/q with mMiss Cut;xB;p/q",numbinxB,binxB,20,0.4,1.4);
      twoD_list.push_back(xB_poq_cutMMiss[i][j]);

      sprintf(temp,"hist_xB_p_%d_%d",i,j);
      hist_xB_p[i][j] = new TH1D(temp,"x_B Distribution;xB;counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_p[i][j]);
    }
  }

  TH1D * hist_xB_q[4];
  for(int i=0; i<4; i++){
      char temp[100];
      sprintf(temp,"hist_xB_q_%d",i);
      hist_xB_q[i] = new TH1D(temp,"x_B Distribution;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_q[i]);
      }


  TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_xB =  new TH1D("hist_xB_Incl" ,"hist;xB;Counts",numbinxBIncl,binxBIncl);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",20,0,5);
  hist_list.push_back(hist_QSq);
  TH2D * xB_QSq = new TH2D("xb_QSq","xB_QSq;xB;QSq;Counts",numbinxB,binxB,40,0,5);
  twoD_list.push_back(xB_QSq);


  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
  }
  for(int i=0; i<twoD_list.size(); i++){
    twoD_list[i]->Sumw2();
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
    
    //Do some things for all events
    if(((i%5000) == 0) && verbose){
      cerr << (i*100.)/(inTree->GetEntries()) <<"% complete \n";
    }
    
    if(nPar > 19){
      cerr<<"There are more than 19 particles in one event! \n Aborting..."<<endl;
      return -1;
    }    
    if(!targInfo.evtxInRange(vtxZCorr[0])){
      continue;
    }
    if(oneSec && (secChoice != getSec(phi))){
	continue;
      }
    


    ////////////////Weight factor
    weight = weight/(targInfo.getLum() * A);

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    //This is for Inclusive measruments only
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    if(TypeData==1){
      //Only look in a range of QSq
      if((QSq<minQSq) || (QSq>=3.5)){ 
	continue;
      }
      if(doMaps){     
	if(targInfo.e_acc(ve) < accMin){
	  continue;
	}
	if(secondTargInfo.e_acc(ve) < accMin){
	  continue;
	}
      weight = weight/targInfo.incl_acc(ve);
      }
      
      hist_xB->Fill(xB,weight);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    //This is for semi-inclusive measruments only
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    else if((TypeData==2) || (TypeData==4)){
      TVector3 vLead = myInfo.getVector(1);
      TVector3 vMiss = vLead - vq;
      double mMiss = myInfo.getMassMiss(1);
      double poq = myInfo.getPoQ(1);
      //double eVTX = eLead - omega + vMiss.Mag();
      //double OS = sq(eVTX) - vMiss.Mag2();
      double thetapq = myInfo.getThetaPQ(1);
      double thetapMq = vq.Angle(vMiss) * (180/M_PI);
      double eff = 1;

      //Get the correct efficiency for the specific particle
      if(myInfo.isProton(1) && doMaps){
	if(targInfo.e_acc(ve)*targInfo.p_acc(vLead) < accMin){
	  continue;
	}
	if(secondTargInfo.e_acc(ve)*secondTargInfo.p_acc(vLead) < accMin){
	  continue;
	}
	//eff = targInfo.e_acc(ve) * targInfo.p_acc(vLead);
	eff = targInfo.semi_acc(ve,vLead);
      }

      weight = weight / eff;

      //Weight related transparency and acceptances
      if(doTrans){
	weight = weight / targInfo.getTrans();
      }


      //Start filling for semi-inclusive
      if( (xB>1.2) && (xB<2) ){
	if( (vMiss.Mag()>0.30) && (vMiss.Mag()<0.60) ){
	  poq_thetapq_0->Fill(poq,thetapq,weight);
	}
	if( (vMiss.Mag()>0.30) && (vMiss.Mag()<1) ){
	  poq_thetapq_1->Fill(poq,thetapq,weight);
	}
      }

      if(thetapq > 25){ continue; }
      if(getPMin(vMiss.Mag()) == -1){ continue; }
      if(getPMax(vMiss.Mag()) == -1){ continue; }
      
      for(int k = getPMin(vMiss.Mag()); k >= 0; k--){
	for(int l = getPMax(vMiss.Mag()); l <= 2; l++){
	  xB_mMiss_All[k][l]->Fill(xB,mMiss,weight);
	  xB_poq_All[k][l]->Fill(xB,poq,weight);
	  if( (poq > 0.62) && (poq < 0.96) ){
	    xB_mMiss_cutPoQ[k][l]->Fill(xB,mMiss,weight);
	  }
	  if( (mMiss > 0.9) && (mMiss < 1.1) ){
	    xB_poq_cutMMiss[k][l]->Fill(xB,poq,weight);
	    hist_xB_p[k][l]->Fill(xB,weight); 
	  }
	}
      }

      if( (mMiss > 0.9) && (mMiss < 1.1) ){
	if( (poq > 0.62) && (poq < 0.96) ){
	  if( (vMiss.Mag() > 0.3) && (vMiss.Mag() < 0.6) ){
	    if( (xB > 1.2) && (xB < 2) ){
	      hist_xB->Fill(xB,weight);	  	  
	    }
	  }
	}
      }

      if( (mMiss > 0.9) && (mMiss < 1.1) ){      
	if( (vMiss.Mag() > 0.3) && (vMiss.Mag() < 0.6) ){
	  
	  xB_QSq->Fill(xB,QSq,weight);    

	  if(QSq > 0){
	    hist_xB_q[0]->Fill(xB,weight);
	  }
	  if(QSq > 1){
	    hist_xB_q[1]->Fill(xB,weight);
	  }
	  if(QSq > 1.25){
	    hist_xB_q[2]->Fill(xB,weight);
	  }
	  if(QSq > 1.5){
	    hist_xB_q[3]->Fill(xB,weight);
	  }
	}
      }

      if(TypeData==4){
	int dIndex = 0;
	for(int j = 2; j < nPar; j++){
	  if(myInfo.isDelta(j)){dIndex=j;}
	}
	if(dIndex != 0){
	  TVector3 vDelta(momx[dIndex],momy[dIndex],momz[dIndex]);
	  TVector3 vCM = (vMiss + vDelta);
	  TVector3 vRel = (vMiss - vDelta);
	
	  //hist_pRel->Fill(vRel.Mag()/2,weight);
	  //hist_pCM->Fill(vCM.Mag(),weight);
	  cout<<vCM.Mag()<<"\n";
	}
      }

      counter++;
    }
    else{
      cerr<<"Your data type("<<TypeData<<") is not in my list.\n Aborting..";
      exit(-2);
    }
    
    //cout<<weight<<"\n";
    //Get cross sections and fill histograms
    hist_Cross->Fill(1,1);
    hist_QSq->Fill(QSq,weight);
    	
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
  for(int i=0; i<twoD_list.size(); i++){
    twoD_list[i]->Write();
  }

  xB_QSq->Write();
  outputFile->Close();
  cerr<< argv[2]<<" has been completed. \n\n\n";
  
  return 0;
}

int getPMin(const double PMiss){

  if(PMiss<0){
    return -1; 
  }
  else if(PMiss<0.1){
    return 0;
  }
  else if(PMiss<0.13){
    return 1;
  }
  else if(PMiss<0.16){
    return 2;
  }
  else if(PMiss<0.2){
    return 3;
  }
  else if(PMiss<0.3){
    return 4;
  }
  else{
    return 5;
  }

}

int getPMax(const double PMiss){

  if(PMiss<0.6){
    return 0;
  }
  else if(PMiss<1){
    return 1;
  }
  else if(PMiss<2){
    return 2;
  }
  else{
    return -1;
  }
}


int getXBBin(const double XB){

  if(XB<1.2){
    return 0;
  }
  else if(XB<1.4){
    return 1;
  }
  else if(XB<1.6){
    return 2;
  }
  else
    return 3;
}

int getPMissBin(const double PMiss){

  if(PMiss<0.15){
    return 0;
  }
  else if(PMiss<0.25){
    return 1;
  }
  else if(PMiss<0.4){
    return 2;
  }
  else
    return 3;
}

int getQSqBin(const double QSq){

  if(QSq<1.7){
    return 0;
  }
  else if(QSq<1.95){
    return 1;
  }
  else if(QSq<2.25){
    return 2;
  }
  else
    return 3;
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

  /*TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_thetaL =  new TH1D("hist_thetaL" ,"hist;thetaL;Counts",180,0,180);
  hist_list.push_back(hist_thetaL);
  TH1D * hist_xB_lowM =  new TH1D("hist_xB_lowM" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB_lowM);
  TH1D * hist_xB_highM =  new TH1D("hist_xB_highM" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB_highM);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",20,0,5);
  hist_list.push_back(hist_QSq);
  TH1D * hist_pMiss_lowM =  new TH1D("hist_pMiss_lowM" ,"hist;pMiss;Counts",50,0,2);
  hist_list.push_back(hist_pMiss_lowM);
  TH1D * hist_pMiss_highM =  new TH1D("hist_pMiss_HighM" ,"hist;pMiss;Counts",50,0,2);
  hist_list.push_back(hist_pMiss_highM);
  TH1D * hist_pRel =  new TH1D("hist_pRel" ,"hist;pRel;Counts",50,0,2);
  hist_list.push_back(hist_pRel);
  TH1D * hist_pCM =  new TH1D("hist_CM" ,"hist;pCM;Counts",50,0,2);
  hist_list.push_back(hist_pCM);
  //Binned histos
  TH1D * hist_xB_binQSq[4];
  TH1D * hist_xB_binPMiss[4];
  TH1D * hist_pMiss_binXB[4];
  for(int i=0; i<4; i++){
      char temp[100];

      sprintf(temp,"hist_xB_bin_QSq%d",i);
      hist_xB_binQSq[i] = new TH1D(temp,"binQSq_hist;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_binQSq[i]);

      sprintf(temp,"hist_xB_bin_PMiss%d",i);
      hist_xB_binPMiss[i] = new TH1D(temp,"binPMiss_hist;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_binPMiss[i]);

      sprintf(temp,"hist_pMiss_bin_XB%d",i);
      hist_pMiss_binXB[i] = new TH1D(temp,"binXB;pMiss;Counts",50,0,2);
      hist_list.push_back(hist_pMiss_binXB[i]);
  }

  //histos with different cuts
  TH1D * hist_xB_cutPMiss[8];
  for(int i=0; i<8; i++){
      char temp[100];
      sprintf(temp,"hist_xB_cut_PMiss%d",i);
      hist_xB_cutPMiss[i] = new TH1D(temp,"cutPMiss_hist;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_cutPMiss[i]);
  }
  TH1D * hist_pMiss_cutXB[4];
  for(int i=0; i<4; i++){
      char temp[100];
      sprintf(temp,"hist_pMiss_cut_XB%d",i);
      hist_pMiss_cutXB[i] = new TH1D(temp,"cutXB_hist;pMiss;Counts",50,0,2);
      hist_list.push_back(hist_pMiss_cutXB[i]);
  }
  TH1D * hist_xB_cutTheta[4];
  for(int i=0; i<4; i++){
      char temp[100];
      sprintf(temp,"hist_xB_cut_Theta_%d",i);
      hist_xB_cutTheta[i] = new TH1D(temp,"cutTheta_hist;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_cutTheta[i]);
      }*/


  /*TH2D * xB_QSq = new TH2D("xb_QSq","xB_QSq;xB;QSq;Counts",numbinxB,binxB,20,0,5);
  twoD_list.push_back(xB_QSq);
  TH2D * xB_pMiss = new TH2D("xb_pMiss","xB_pMiss;xB;pMiss;Counts",numbinxB,binxB,50,0,2);
  twoD_list.push_back(xB_pMiss);
  TH2D * eVTX_pMiss = new TH2D("eVTX_pMiss","eVTX_pMiss;eVTX;pMiss;Counts",40,-0.5,3.5,50,0,2);
  twoD_list.push_back(eVTX_pMiss);
  TH2D * OS_pMiss = new TH2D("OS_pMiss","OS_pMiss;OS;pMiss;Counts",40,-0.5,1.5,50,0,2);
  twoD_list.push_back(OS_pMiss);
  TH2D * pMiss_poq = new TH2D("pMiss_poq","pMiss_pq;pMiss;poq;Counts",50,0,2,28,0,1.4);
  twoD_list.push_back(pMiss_poq);
  TH2D * xB_thetaL = new TH2D("xb_thetaL","xB_thetaL;xB;thetaL;Counts",numbinxB,binxB,180,0,180);
  twoD_list.push_back(xB_thetaL);
  TH2D * thetaL_mMiss = new TH2D("xb_mMiss","thetaL_mMiss;thetaL;mMiss;Counts",180,0,180,18,0.6,1.5);
  twoD_list.push_back(thetaL_mMiss);
  TH2D * xB_thetapq = new TH2D("xB_thetapq","xB_pq;xB;thetapq;Counts",numbinxB,binxB,35,0,35);
  twoD_list.push_back(xB_thetapq);
  TH2D * xB_thetapMq = new TH2D("xB_thetapMq","xB_pMq;xB;thetapMq;Counts",numbinxB,binxB,180,0,180);
  twoD_list.push_back(xB_thetapMq);
  TH2D * thetapMq_thetapq = new TH2D("thetapMq_thetapq","thetapMq_pq;thetapMq;thetapq;Counts",180,0,180,35,0,35);
  twoD_list.push_back(thetapMq_thetapq);
  TH2D * pMiss_thetapq = new TH2D("pMiss_thetapq","pMiss_pq;pMiss;thetapq;Counts",50,0,2,35,0,35);
  twoD_list.push_back(pMiss_thetapq);
  TH2D * pLead_thetapq = new TH2D("pLead_thetapq","pLead_pq;pLead;thetapq;Counts",20,0,2,35,0,35);
  twoD_list.push_back(pLead_thetapq);
  TH2D * pMT_thetapq = new TH2D("pMT_thetapq","pMT_pq;pMT;thetapq;Counts",40,-2,2,35,0,35);
  twoD_list.push_back(pMT_thetapq);
  TH2D * pLT_thetapq = new TH2D("pLT_thetapq","pLT_pq;pLT;thetapq;Counts",40,0,4,35,0,35);
  twoD_list.push_back(pLT_thetapq);
  TH2D * plead_poq_cutPMiss[8];
  TH2D * xB_eVTX_cutPMiss[8];
  TH2D * xB_OS_cutPMiss[8];
  for(int i=0; i<8; i++){
      char temp[100];

      sprintf(temp,"xB_mMiss_cut_PMiss%d",i);
      xB_mMiss_cutPMiss[i] = new TH2D(temp,"x_B vs m_Miss;xB;mMiss",numbinxB,binxB,50,0.5,3.5);
      twoD_list.push_back(xB_mMiss_cutPMiss[i]);

      sprintf(temp,"xB_poq_cut_PMiss%d",i);
      xB_poq_cutPMiss[i] = new TH2D(temp,"x_B vs p/q;xB;p/q",numbinxB,binxB,28,0,1.4);
      twoD_list.push_back(xB_poq_cutPMiss[i]);

      sprintf(temp,"plead_poq_cut_PMiss%d",i);
      plead_poq_cutPMiss[i] = new TH2D(temp,"p_lead vs p/q;plead;p/q",100,0,4,28,0,1.4);
      twoD_list.push_back(plead_poq_cutPMiss[i]);

      sprintf(temp,"xB_eVTX_cut_PMiss%d",i);
      xB_eVTX_cutPMiss[i] = new TH2D(temp,"x_B vs eVTX;xB;eVTX",numbinxB,binxB,40,-0.5,3.5);
      twoD_list.push_back(xB_eVTX_cutPMiss[i]);

      sprintf(temp,"xB_OS_cut_PMiss%d",i);
      xB_OS_cutPMiss[i] = new TH2D(temp,"x_B vs OS;xB;OS",numbinxB,binxB,40,-0.5,1.5);
      twoD_list.push_back(xB_OS_cutPMiss[i]);

  }
  */


      //hist_thetaL->Fill((180/M_PI)*vLead.Theta(),weight);
      //xB_thetaL->Fill(xB,(180/M_PI)*vLead.Theta(),weight);
      //thetaL_mMiss->Fill((180/M_PI)*vLead.Theta(),mMiss,weight);

      //if( (vMiss.Mag()>2) || (vMiss.Mag()<0.3) ){ continue; }
      //if(QSq < 1.5){ continue; }
      //if(eVTX > 1.75){ continue; }
      //if((mMiss < 0.9) || (mMiss > 1.1)){ continue; }

      /*
      pMiss_poq->Fill(vMiss.Mag(),poq,weight);
      xB_pMiss->Fill(xB,vMiss.Mag(),weight);
      eVTX_pMiss->Fill(eVTX,vMiss.Mag(),weight);
      OS_pMiss->Fill(OS,vMiss.Mag(),weight);
      hist_xB_binPMiss[getPMissBin(vMiss.Mag())]->Fill(xB,weight);
      hist_pMiss_binXB[getXBBin(xB)]->Fill(vMiss.Mag(),weight);
      
      if(vMiss.Mag()>0.00){
	hist_xB_cutPMiss[0]->Fill(xB,weight);
	xB_mMiss_cutPMiss[0]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[0]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[0]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[0]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[0]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.10){
	hist_xB_cutPMiss[1]->Fill(xB,weight);
	xB_mMiss_cutPMiss[1]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[1]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[1]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[1]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[1]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.13){
	hist_xB_cutPMiss[2]->Fill(xB,weight);
	xB_mMiss_cutPMiss[2]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[2]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[2]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[2]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[2]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.16){
	hist_xB_cutPMiss[3]->Fill(xB,weight);
	xB_mMiss_cutPMiss[3]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[3]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[3]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[3]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[3]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.20){
	hist_xB_cutPMiss[4]->Fill(xB,weight);
	xB_mMiss_cutPMiss[4]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[4]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[4]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[4]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[4]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.30){
	hist_xB_cutPMiss[5]->Fill(xB,weight);
	xB_mMiss_cutPMiss[5]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[5]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[5]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[5]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[5]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.40){
	hist_xB_cutPMiss[6]->Fill(xB,weight);
	xB_mMiss_cutPMiss[6]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[6]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[6]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[6]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[6]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.50){
	hist_xB_cutPMiss[7]->Fill(xB,weight);
	xB_mMiss_cutPMiss[7]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[7]->Fill(xB,poq,weight);
	plead_poq_cutPMiss[7]->Fill(vLead.Mag(),poq,weight);
	xB_eVTX_cutPMiss[7]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[7]->Fill(xB,OS,weight);
      }
     

      if(xB>0.8){hist_pMiss_cutXB[0]->Fill(vMiss.Mag());}
      if(xB>1.0){hist_pMiss_cutXB[1]->Fill(vMiss.Mag());}
      if(xB>1.2){hist_pMiss_cutXB[2]->Fill(vMiss.Mag());}
      if(xB>1.4){hist_pMiss_cutXB[3]->Fill(vMiss.Mag());}
    
      if(vMiss.Mag()>0.3){
	xB_QSq->Fill(xB,QSq,weight);    
	hist_xB_binQSq[getQSqBin(QSq)]->Fill(xB,weight);
      }

      if( vMiss.Mag() < 0.3){ continue; }
      if( thetapMq < 40){ continue; }

      if((mMiss > 0.9)&&(mMiss < 1.1)){
	hist_xB_lowM->Fill(xB,weight);
	hist_pMiss_lowM->Fill(vMiss.Mag(),weight);      
      }
      else if((mMiss > 1.2)&&(mMiss < 1.45)){
	hist_xB_highM->Fill(xB,weight);
	hist_pMiss_highM->Fill(vMiss.Mag(),weight);      
      }


      xB_thetapq->Fill(xB,thetapq,weight);
      xB_thetapMq->Fill(xB,thetapMq,weight);
      thetapMq_thetapq->Fill(thetapMq,thetapq,weight);
      pMiss_thetapq->Fill(vMiss.Mag(),thetapq,weight);
      pLead_thetapq->Fill(vLead.Mag(),thetapq,weight);
      pMT_thetapq->Fill(vMiss.Dot(vq.Unit()),thetapq,weight);
      pLT_thetapq->Fill(vLead.Dot(vq.Unit()),thetapq,weight);

      if(thetapq < 25){hist_xB_cutTheta[0]->Fill(xB,weight);}
      if(thetapq < 15){hist_xB_cutTheta[1]->Fill(xB,weight);}
      if(thetapq < 11){hist_xB_cutTheta[2]->Fill(xB,weight);}
      if(thetapq < 7){hist_xB_cutTheta[3]->Fill(xB,weight);}
      */
