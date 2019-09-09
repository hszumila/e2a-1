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
const TVector3 vBeam(0.,0.,4.461);
double sq(double x){
  return x*x;
}
int getQSqBin(const double QSq);
int getXBBin(const double XB);
int getPMissBin(const double PMiss);
int getSec(const double phi);

void help_message()
{
  cerr<< "Argumets: ./get_hist /path/to/input/skim/file /path/to/output/file [Nucleus A] [Data Type] [optional flags]\n\n"
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
  int TypeData = atoi(argv[4]);
  //Get target info
  target_Info targInfo(A);  
  
  bool verbose = false;
  bool doMaps = true;
  bool doTrans = true;
  bool oneSec = false;
  double minQSq = 1.4;
  int secChoice = -1;

  
  int c;
  while ((c=getopt (argc-4, &argv[4], "hvn:x:ms:tQ:")) != -1) //First two arguments are not optional flags.
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
  double binxB[] = {0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.38,1.58,1.79,2};
  int numbinxB = 22;
  //double binxB[] = {0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.26,1.32,1.38,1.44,1.53,1.62,1.72,1.8,1.9,2};
  //int numbinxB = 26;

  //Make a histogram list to make things easier
  vector<TH2*> twoD_list;
  TH2D * xB_QSq = new TH2D("xb_QSq","xB_QSq;xB;QSq;Counts",numbinxB,binxB,20,0,5);
  twoD_list.push_back(xB_QSq);
  TH2D * xB_pMiss = new TH2D("xb_pMiss","xB_pMiss;xB;pMiss;Counts",numbinxB,binxB,15,0,0.6);
  twoD_list.push_back(xB_pMiss);
  TH2D * eVTX_pMiss = new TH2D("eVTX_pMiss","eVTX_pMiss;eVTX;pMiss;Counts",40,-0.5,1.5,15,0,0.6);
  twoD_list.push_back(eVTX_pMiss);
  TH2D * OS_pMiss = new TH2D("OS_pMiss","OS_pMiss;OS;pMiss;Counts",40,-0.5,1.5,15,0,0.6);
  twoD_list.push_back(OS_pMiss);
  TH2D * xB_mMiss_cutPMiss[6];
  TH2D * xB_poq_cutPMiss[6];
  TH2D * xB_eVTX_cutPMiss[6];
  TH2D * xB_OS_cutPMiss[6];
  for(int i=0; i<6; i++){
      char temp[100];

      sprintf(temp,"xB_mMiss_cut_PMiss%d",i);
      xB_mMiss_cutPMiss[i] = new TH2D(temp,"x_B vs m_Miss;xB;mMiss",numbinxB,binxB,18,0.6,1.5);
      twoD_list.push_back(xB_mMiss_cutPMiss[i]);

      sprintf(temp,"xB_poq_cut_PMiss%d",i);
      xB_poq_cutPMiss[i] = new TH2D(temp,"x_B vs p/q;xB;p/q",numbinxB,binxB,28,0,1.4);
      twoD_list.push_back(xB_poq_cutPMiss[i]);

      sprintf(temp,"xB_eVTX_cut_PMiss%d",i);
      xB_eVTX_cutPMiss[i] = new TH2D(temp,"x_B vs eVTX;xB;eVTX",numbinxB,binxB,40,-0.5,1.5);
      twoD_list.push_back(xB_eVTX_cutPMiss[i]);

      sprintf(temp,"xB_OS_cut_PMiss%d",i);
      xB_OS_cutPMiss[i] = new TH2D(temp,"x_B vs OS;xB;OS",numbinxB,binxB,40,-0.5,1.5);
      twoD_list.push_back(xB_OS_cutPMiss[i]);

  }
  for(int i=0; i<twoD_list.size(); i++){
    twoD_list[i]->Sumw2();
  }


  vector<TH1*> hist_list;
  TH1D * hist_Cross = new TH1D("totalCross","Cross;1;Value",1,0.5,1.5);
  hist_list.push_back(hist_Cross);
  TH1D * hist_xB_lowM =  new TH1D("hist_xB_lowM" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB_lowM);
  TH1D * hist_xB_highM =  new TH1D("hist_xB_highM" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB_highM);
  TH1D * hist_xB =  new TH1D("hist_xB" ,"hist;xB;Counts",numbinxB,binxB);
  hist_list.push_back(hist_xB);
  TH1D * hist_QSq =  new TH1D("hist_QSq" ,"hist;QSq;Counts",20,0,5);
  hist_list.push_back(hist_QSq);
  TH1D * hist_pMiss_lowM =  new TH1D("hist_pMiss_lowM" ,"hist;pMiss;Counts",15,0,0.6);
  hist_list.push_back(hist_pMiss_lowM);
  TH1D * hist_pMiss_highM =  new TH1D("hist_pMiss_HighM" ,"hist;pMiss;Counts",15,0,0.6);
  hist_list.push_back(hist_pMiss_highM);
  TH1D * hist_pRel =  new TH1D("hist_pRel" ,"hist;pRel;Counts",15,0,0.6);
  hist_list.push_back(hist_pRel);
  TH1D * hist_pCM =  new TH1D("hist_CM" ,"hist;pCM;Counts",15,0,0.6);
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
      hist_pMiss_binXB[i] = new TH1D(temp,"binXB;pMiss;Counts",15,0,0.6);
      hist_list.push_back(hist_pMiss_binXB[i]);
  }

  //histos with different cuts
  TH1D * hist_xB_cutPMiss[6];
  for(int i=0; i<6; i++){
      char temp[100];
      sprintf(temp,"hist_xB_cut_PMiss%d",i);
      hist_xB_cutPMiss[i] = new TH1D(temp,"cutPMiss_hist;xB;Counts",numbinxB,binxB);
      hist_list.push_back(hist_xB_cutPMiss[i]);
  }
  TH1D * hist_pMiss_cutXB[4];
  for(int i=0; i<4; i++){
      char temp[100];
      sprintf(temp,"hist_pMiss_cut_XB%d",i);
      hist_pMiss_cutXB[i] = new TH1D(temp,"cutXB_hist;pMiss;Counts",15,0,0.6);
      hist_list.push_back(hist_pMiss_cutXB[i]);
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


  //Loop over TTree
  for(int i = 0; i < inTree->GetEntries(); i++){
    double weight = 1;
    inTree->GetEntry(i);
    event_Info myInfo(nPar,parID,xB,momx,momy,momz,vtxZCorr,0,0,false,false);

    TVector3 ve(momx[0],momy[0],momz[0]);
    TVector3 vq = vBeam - ve;
    double omega = vBeam.Mag() - ve.Mag();
    double phi = ve.Phi() * 180 / M_PI;
    //hist_eAccu->Fill(targInfo.e_acc(ve));
    
    //Display completed
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
    weight = 1/(targInfo.getLum() * A);
    if(doMaps){     
      if(targInfo.e_acc(ve) < 0.4){
	continue;
      }
      weight = weight/targInfo.e_acc(ve);
    }
  
    //Now do some calculations only if you have the correct number of particles
    if(TypeData==1){
      //Only look in a range of QSq
      if((QSq<minQSq) || (QSq>=3.5)){ 
	continue;
      }
    }
    else if((TypeData==2) || (TypeData==4)){
      TVector3 vLead(momx[1],momy[1],momz[1]);
      TVector3 vMiss = vLead - vq;
      double eLead = sqrt(vLead.Mag2() + sq(mP));
      double eMiss = omega + mP + mN - eLead;
      double mMiss = sqrt(sq(eMiss)-vMiss.Mag2());
      double poq = vLead.Mag()/vq.Mag();
      double eVTX = eLead - omega;
      double OS = sq(eVTX) - vMiss.Mag2();
      //sqrt(vMiss.Mag2() - sq(mP))      
      //      cout<<eVTX<<"\n";
      double eff = 0;
      //Get the correct efficiency for the specific particle
      if(myInfo.isProton(1)){
	eff = targInfo.p_acc(vLead);
      }
      else if(myInfo.isNeutron(1)){
	eff = 1;
      }
      
      //Weight related transparency and acceptances
      if(doTrans){
	weight = weight / targInfo.getTrans();
      }
      if(doMaps){
	if(targInfo.p_acc(vLead) < 0.4){
	  continue;
	}      
	weight = weight / eff;
      }
      
      if((mMiss > 0.9)&&(mMiss < 1.1)){
	hist_xB_lowM->Fill(xB,weight);
	hist_pMiss_lowM->Fill(vMiss.Mag(),weight);      
      }
      else if((mMiss > 1.2)&&(mMiss < 1.45)){
	hist_xB_highM->Fill(xB,weight);
	hist_pMiss_highM->Fill(vMiss.Mag(),weight);      
      }
      
      if((mMiss < 0.9) || (mMiss > 1.1)){ continue; }
      //if((eVTX < 0)){ continue; }
    
      xB_pMiss->Fill(xB,vMiss.Mag(),weight);
      eVTX_pMiss->Fill(eVTX,vMiss.Mag(),weight);
      OS_pMiss->Fill(OS,vMiss.Mag(),weight);
      hist_xB_binPMiss[getPMissBin(vMiss.Mag())]->Fill(xB,weight);
      hist_pMiss_binXB[getXBBin(xB)]->Fill(vMiss.Mag(),weight);

      
      if(vMiss.Mag()>0.00){
	hist_xB_cutPMiss[0]->Fill(xB,weight);
	xB_mMiss_cutPMiss[0]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[0]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[0]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[0]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.10){
	hist_xB_cutPMiss[1]->Fill(xB,weight);
	xB_mMiss_cutPMiss[1]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[1]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[1]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[1]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.13){
	hist_xB_cutPMiss[2]->Fill(xB,weight);
	xB_mMiss_cutPMiss[2]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[2]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[2]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[2]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.16){
	hist_xB_cutPMiss[3]->Fill(xB,weight);
	xB_mMiss_cutPMiss[3]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[3]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[3]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[3]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.20){
	hist_xB_cutPMiss[4]->Fill(xB,weight);
	xB_mMiss_cutPMiss[4]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[4]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[4]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[4]->Fill(xB,OS,weight);
      }
      if(vMiss.Mag()>0.30){
	hist_xB_cutPMiss[5]->Fill(xB,weight);
	xB_mMiss_cutPMiss[5]->Fill(xB,mMiss,weight);
	xB_poq_cutPMiss[5]->Fill(xB,poq,weight);
	xB_eVTX_cutPMiss[5]->Fill(xB,eVTX,weight);
	xB_OS_cutPMiss[5]->Fill(xB,OS,weight);
      }
     

      if(xB>0.8){hist_pMiss_cutXB[0]->Fill(vMiss.Mag());}
      if(xB>1.0){hist_pMiss_cutXB[1]->Fill(vMiss.Mag());}
      if(xB>1.2){hist_pMiss_cutXB[2]->Fill(vMiss.Mag());}
      if(xB>1.4){hist_pMiss_cutXB[3]->Fill(vMiss.Mag());}
    
      if(vMiss.Mag()>0.3){
	xB_QSq->Fill(xB,QSq,weight);    
	hist_xB_binQSq[getQSqBin(QSq)]->Fill(xB,weight);
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
	
	  hist_pRel->Fill(vRel.Mag()/2,weight);
	  hist_pCM->Fill(vCM.Mag(),weight);
	  cout<<vCM.Mag()<<"\n";
	}
      }
    }
    else{
      cerr<<"Your data type("<<TypeData<<") is not in my list.\n Aborting..";
      exit(-2);
    }

    //Get cross sections and fill histograms
    hist_Cross->Fill(1,weight);
    hist_xB->Fill(xB,weight);
    hist_QSq->Fill(QSq,weight);
    xB_QSq->Fill(xB,QSq,weight);    
	
  }
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
