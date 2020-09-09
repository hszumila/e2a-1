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
#include "TGraph.h"
#include "TVector3.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "event_Info.h"
#include "target_Info.h"
#include "e2a_constants.h"

using namespace std;

double sq(double x){ return x*x; };
int getSec(const double phi);

void help_message()
{
  cerr<< "Argumets: ./vtxCheck /path/to/input/skim/file /path/to/output/root/file /output/path/to/electron/angle/text/file /output/path/to/proton/angle/text/file A\n\n\n";
}

int main(int argc, char ** argv){
     
  if (argc != 6)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_message();
      return -1;
    } 
  
  //Read in arguments.
  TFile * inputFile = new TFile(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");
  ofstream electrontextFile;
  electrontextFile.open(argv[3]);
  ofstream protontextFile;
  protontextFile.open(argv[4]);
  int A = atoi(argv[5]);

  TH1D * hist_eVTX =  new TH1D("hist_eVTX" ,"hist;eVTX;Counts",400,-10,10);


  vector<TH2*> hist_list;
  TH2D * hist_eVTX_pVTX =  new TH2D("hist_eVTX_pVTX" ,"hist;eVTX;pVTX;Counts",400,-10,10,400,-10,10);
  hist_list.push_back(hist_eVTX_pVTX);
  TH2D * hist_eVTX_eTheta =  new TH2D("hist_eVTX_eTheta" ,"hist;eVTX;eTheta;Counts",400,-10,10,180,0,180);
  hist_list.push_back(hist_eVTX_eTheta);
  TH2D * hist_pVTX_pTheta =  new TH2D("hist_pVTX_pTheta" ,"hist;pVTX;pTheta;Counts",400,-10,10,180,0,180);
  hist_list.push_back(hist_pVTX_pTheta);
  TH2D * hist_eVTX_eMom =  new TH2D("hist_eVTX_eMom" ,"hist;eVTX;eMom;Counts",400,-10,10,100,0,6);
  hist_list.push_back(hist_eVTX_eMom);
  TH2D * hist_pVTX_pMom =  new TH2D("hist_pVTX_pMom" ,"hist;pVTX;pMom;Counts",400,-10,10,100,0,6);
  hist_list.push_back(hist_pVTX_pMom);
  TH2D * hist_diffVTX_eTheta =  new TH2D("hist_diffVTX_eTheta" ,"hist;diffVTX;eTheta;Counts",400,-10,10,180,0,180);
  hist_list.push_back(hist_diffVTX_eTheta);
  TH2D * hist_diffVTX_pTheta =  new TH2D("hist_diffVTX_pTheta" ,"hist;diffVTX;pTheta;Counts",400,-10,10,180,0,180);
  hist_list.push_back(hist_diffVTX_pTheta);
  TH2D * hist_diffVTX_eMom =  new TH2D("hist_diffVTX_eMom" ,"hist;diffVTX;eMom;Counts",400,-10,10,100,0,6);
  hist_list.push_back(hist_diffVTX_eMom);
  TH2D * hist_diffVTX_pMom =  new TH2D("hist_diffVTX_pMom" ,"hist;diffVTX;pMom;Counts",400,-10,10,100,0,6);
  hist_list.push_back(hist_diffVTX_pMom);


  TH2D * hist_eVTX_pVTX_sec[6];
  TH2D * hist_diffVTX_eTheta_sec[6];
  TH2D * hist_diffVTX_pTheta_sec[6];
  for(int i=0; i<6; i++){
    char temp[100];
    sprintf(temp,"hist_eVTX_pVTX_sec_%d",i);
    hist_eVTX_pVTX_sec[i] =  new TH2D(temp,"hist;eVTX;pVTX;Counts",400,-10,10,400,-10,10);
    hist_list.push_back(hist_eVTX_pVTX_sec[i]);

    sprintf(temp,"hist_diffVTX_eTheta_sec_%d",i);
    hist_diffVTX_eTheta_sec[i] =  new TH2D(temp ,"hist;diffVTX;eTheta;Counts",400,-10,10,164,16,180);
    hist_list.push_back(hist_diffVTX_eTheta_sec[i]);

    sprintf(temp,"hist_diffVTX_pTheta_sec_%d",i);
    hist_diffVTX_pTheta_sec[i] =  new TH2D(temp ,"hist;diffVTX;pTheta;Counts",400,-10,10,164,16,180);
    hist_list.push_back(hist_diffVTX_pTheta_sec[i]);
  }


  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
  }

  target_Info targInfo(A);  

  //Make trees and histograms for data
  TTree * inTree = (TTree*)inputFile->Get("T");

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


  cerr<<"input TTree opened from input file \n output TTree creadted. \n";
  
  int fin = inTree->GetEntries();
  int nEvents = 0;

  //Loop over TTree
  for(int i = 0; i < fin; i++){

    inTree->GetEntry(i);
    
    //Display completed
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }    
    
    if(nPar > 2){
      continue;
    }

    event_Info myInfo(nPar,parID,xB,QSq,momx,momy,momz,vtxZCorr);
    double eMom = myInfo.getVector(0).Mag();
    double pMom = myInfo.getVector(1).Mag(); 
    double ePhi = myInfo.getVector(0).Phi() * 180 / M_PI;
    double eTheta = myInfo.getVector(0).Theta() * 180 / M_PI;
    double pTheta = myInfo.getVector(1).Theta() * 180 / M_PI; 
    double eVTX = vtxZCorr[0];
    double pVTX = vtxZCorr[1];
    double diffVTX = eVTX - pVTX;

    //if(!targInfo.vtxInRange(eVTX,pVTX,eTheta,pTheta)){
    //  continue;
    //}

    hist_eVTX->Fill(eVTX);
    hist_eVTX_pVTX->Fill(eVTX,pVTX);
    hist_eVTX_pVTX_sec[getSec(ePhi)]->Fill(eVTX,pVTX);



    if(!targInfo.semiVTXInRange(myInfo,1)){
      continue;
    }

    hist_diffVTX_eTheta_sec[getSec(ePhi)]-> Fill(diffVTX ,eTheta);
    hist_diffVTX_pTheta_sec[getSec(ePhi)]-> Fill(diffVTX ,pTheta);

    hist_eVTX_eTheta->    Fill(eVTX    ,eTheta);
    hist_pVTX_pTheta->    Fill(pVTX    ,pTheta);
    hist_eVTX_eMom->      Fill(eVTX    ,eMom  );
    hist_pVTX_pMom->      Fill(pVTX    ,pMom  );
    hist_diffVTX_eTheta-> Fill(diffVTX ,eTheta);
    hist_diffVTX_pTheta-> Fill(diffVTX ,pTheta);
    hist_diffVTX_eMom  -> Fill(diffVTX ,eMom  );
    hist_diffVTX_pMom  -> Fill(diffVTX ,pMom  );

    nEvents++;
  }

  /*
  for(int l = 0; l < 6; l++){
    char temp[100];

    sprintf(temp,"eGraph_Left_sec_%d",l);
    eGraphLeft[k] = new TGraphAsymmErrors();
    eGraphLeft[k]->SetName(temp);

    sprintf(temp,"eGraph_Right_sec_%d",l);
    eGraphRight[k] = new TGraphAsymmErrors();
    eGraphRight[k]->SetName(temp);

    sprintf(temp,"pGraphp_Left_sec_%d",l);
    pGraphLeft[k] = new TGraphAsymmErrors();
    pGraphLeft[k]->SetName(temp);

    sprintf(temp,"pGraph_Right_sec_%d",l);
    pGraphRight[k] = new TGraphAsymmErrors();
    pGraphRight[k]->SetName(temp);
  }  

  
  for(int k = 0; k < 164; k+=2){
    double eY1 = hist_diffVTX_eTheta->GetYaxis()->GetBinCenter(k);
    double eY2 = hist_diffVTX_eTheta->GetYaxis()->GetBinCenter(k+1);
    double eY = (eY1+eY2)/2;
    electrontextFile<<eY<<" ";
    for(int j = 0; j<6; j++){
      TH1D * e_proj = hist_diffVTX_eTheta_sec[j]->ProjectionX("sliceProj",k,k+2);
      double eNum= e_proj->GetEntries();
      double eMax = e_proj->GetMaximum();
      double eMean = e_proj->GetMean();
      double eSTD = e_proj->GetStdDev();
      if(eNum>100){
	TF1 * eFit = new TF1("myFit","gaus",-5,5);
	eFit->SetParameter(0,eMax);
	eFit->SetParameter(1,eMean);
	eFit->SetParameter(2,eSTD);
	TFitResultPtr ePoint = e_proj->Fit(eFit,"qesrn","",-5,5);
	double eXL=0;
	double eXR=0;
	if(ePoint == 0){
	  eMax = ePoint->Parameter(0);
	  eMean = ePoint->Parameter(1);
	  eSTD = ePoint->Parameter(2);
	  eXL = eMean - (2 * eSTD);
	  eXR = eMean + (2 * eSTD);
	  eGraphLeft->SetPoint(eGraphLeft->GetN(),eXL,eY);
	  eGraphRight->SetPoint(eGraphRight->GetN(),eXR,eY);
	  electrontextFile<<eXL<<" "<<eXR<<" ";
	}	
      }
    }
  }
  */

  int eN = hist_diffVTX_eTheta->GetYaxis()->GetNbins();
  TGraph * eGraphLeft = new TGraph();
  eGraphLeft->SetName("eGraph_Left");
  TGraph * eGraphRight = new TGraph();
  eGraphRight->SetName("eGraph_Right");
  for(int k = 0; k < eN; k+=2){
    TH1D * e_proj = hist_diffVTX_eTheta->ProjectionX("sliceProj",k,k+2);
    double eNum= e_proj->GetEntries();
    double eMax = e_proj->GetMaximum();
    double eMean = e_proj->GetMean();
    double eSTD = e_proj->GetStdDev();
    double eY1 = hist_diffVTX_eTheta->GetYaxis()->GetBinCenter(k);
    double eY2 = hist_diffVTX_eTheta->GetYaxis()->GetBinCenter(k+1);
    double eY = (eY1+eY2)/2;
    if(eNum>100){
      TF1 * eFit = new TF1("myFit","gaus",-5,5);
      eFit->SetParameter(1,eMean);
      eFit->SetParameter(2,eSTD);
      TFitResultPtr ePoint = e_proj->Fit(eFit,"qesrn","",-5,5);
      double eXL=0;
      double eXR=0;
      if(ePoint == 0){
	eMax = ePoint->Parameter(0);
	eMean = ePoint->Parameter(1);
	eSTD = ePoint->Parameter(2);
	eXL = eMean - (2 * eSTD);
	eXR = eMean + (2 * eSTD);
	eGraphLeft->SetPoint(eGraphLeft->GetN(),eXL,eY);
	eGraphRight->SetPoint(eGraphRight->GetN(),eXR,eY);
	electrontextFile<<eY<<" "<<eXL<<" "<<eXR<<" \n";
      }
    }
  }

  int pN = hist_diffVTX_pTheta->GetYaxis()->GetNbins();
  TGraph * pGraphLeft = new TGraph();
  pGraphLeft->SetName("pGraph_Left");
  TGraph * pGraphRight = new TGraph();
  pGraphRight->SetName("pGraph_Right");
  for(int k = 0; k < pN; k+=5){
    TH1D * p_proj = hist_diffVTX_pTheta->ProjectionX("sliceProj",k,k+5);
    double pNum= p_proj->GetEntries();
    double pMax = p_proj->GetMaximum();
    double pMean = p_proj->GetMean();
    double pSTD = p_proj->GetStdDev();
    double pY1 = hist_diffVTX_pTheta->GetYaxis()->GetBinCenter(k);
    double pY2 = hist_diffVTX_pTheta->GetYaxis()->GetBinCenter(k+4);
    double pY = (pY1+pY2)/2;
    if(pNum>100){
      TF1 * pFit = new TF1("myFit","gaus",-5,5);
      pFit->SetParameter(1,pMean);
      pFit->SetParameter(2,pSTD);
      TFitResultPtr pPoint = p_proj->Fit(pFit,"qesrn","",-5,5);
      double pXL=0;
      double pXR=0;
      if(pPoint == 0){
	pMax = pPoint->Parameter(0);
	pMean = pPoint->Parameter(1);
	pSTD = pPoint->Parameter(2);
	pXL = pMean - (2 * pSTD);
	pXR = pMean + (2 * pSTD);
	pGraphLeft->SetPoint(pGraphLeft->GetN(),pXL,pY);
	pGraphRight->SetPoint(pGraphRight->GetN(),pXR,pY);
	protontextFile<<pY<<" "<<pXL<<" "<<pXR<<" \n";
      }
    }
  }


  cout<<"Finished filling tree for output"<<endl;
  cout<<"Number of Events with lead Nucleons: "<<nEvents<<endl;
  outputFile->cd();

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  hist_eVTX->Write();
  eGraphLeft->Write();
  eGraphRight->Write();
  pGraphLeft->Write();
  pGraphRight->Write();
  inputFile->Close();
  outputFile->Close();
  electrontextFile.close();
  protontextFile.close();
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


	/*double root = -2 * log(eMinValue/eMax);
	if(root > 0){
	  eXL = eMean - (eSTD * sqrt(root));
	  eXR = eMean + (eSTD * sqrt(root));
	}
	if(eXL > -0.5){eXL = -0.5;}
	if(eXR <  0.5){eXR =  0.5;}*/

	/*double root = -2 * log(pMinValue/pMax);
	if(root > 0){
	  pXL = pMean - (pSTD * sqrt(root));
	  pXR = pMean + (pSTD * sqrt(root));
	}
	if(pXL > -0.5){pXL = -0.5;}
	if(pXR <  0.5){pXR =  0.5;}*/
