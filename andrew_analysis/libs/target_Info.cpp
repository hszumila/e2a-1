#include "target_Info.h"
#include "Acceptance.h"
target_Info::target_Info(int A)
{
  e2adir = std::string(getenv("E2A_INSTALL"));  
  M_a=(double)A;
  
  if(A==3){
    massA = 2.808392;
    trans = 0.82;
    vzMax = -0.2;
    vzMin = -2.8;
    density = 0.0655;
    ltcc=10224579.8;
    thick=vzMax-vzMin;
    acc_Name = "He3";
    fid_Name = "3He";
    //No 3He rad correction
    rad_Name = "4He";
    vtx_Name = fid_Name;
  }
  else if(A==4){
    massA = 3.727379;
    trans = 0.75;
    vzMax = 1.5;
    vzMin = -2.4;
    density = 0.1370;
    ltcc=6388866.74;
    thick=vzMax-vzMin;
    acc_Name = "He4";
    fid_Name = "4He";
    rad_Name = fid_Name;
    vtx_Name = fid_Name;
  }
  else if(A==12){
    massA = 11.1748632;
    trans = 0.53;
    vzMax = 6;
    vzMin = 4;
    //    vzMax = 5.75;
    //    vzMin = 4.25;
    density = 1.786;
    ltcc=17962592.69;
    thick=0.1;
    acc_Name = "solid";
    fid_Name = "12C";
    rad_Name = fid_Name;
    vtx_Name = fid_Name;
  }
  else{
    std::cerr << "The nucleus you have chosen could not be found\n\n Aborting...\n\n";
    exit(-2);
  }
  eMap = new Acceptance(acc_Name,4461,2250,"e");
  pMap = new Acceptance(acc_Name,4461,2250,"p");
  pipMap = new Acceptance(acc_Name,4461,2250,"pip");
  targFid = new Fiducial(4461,2250,5996,fid_Name,true);
  fillRadArray();
  fillVTXArrayElectron();
  fillVTXArrayProton();
  setLum();

}

target_Info::~target_Info()
{
}

double target_Info::incl_acc(const TVector3 ve)
{
  double sumAcc = 0;
  TRandom3 myRand(0);
  for(int i = 0; i < 100; i++){
    double phiE = 2 * M_PI * myRand.Rndm();
    TVector3 vePrime = ve;
    vePrime.Rotate(phiE,vBeam);
    if(pass_incl_fid(vePrime)){
      sumAcc += e_acc(vePrime);
    } 
  }
  return sumAcc/100;
}

double target_Info::semi_acc(const TVector3 ve,const TVector3 vLead)
{
  double sumAcc = 0;
  TRandom3 myRand(0);

  for(int i = 0; i < 1000; i++){
    double phiE = 2 * M_PI * myRand.Rndm();
    double phiP = 2 * M_PI * myRand.Rndm();

    TVector3 vePrime = ve;
    TVector3 vLeadPrime = vLead;
    vePrime.Rotate(phiE,vBeam);
    vLeadPrime.Rotate(phiE,vBeam);
    TVector3 q = vBeam - vePrime;
    vLeadPrime.Rotate(phiP,q);

    //Before adding it to the average, check that it passed the fiducial
    if(pass_semi_fid(vePrime,vLeadPrime)){
      sumAcc += (e_acc(vePrime) * p_acc(vLeadPrime));
    }
  }
  return sumAcc/1000;

}

bool target_Info::pass_incl_fid(const TVector3 ve)
{
  return targFid->e_inFidRegion(ve);
}

bool target_Info::pass_semi_fid(const TVector3 ve, const TVector3 vLead)
{
  if(targFid->e_inFidRegion(ve) && targFid->pFiducialCut(vLead)){
    return true;
  }
  return false; 
}

double target_Info::e_acc(TVector3 p)
{
  return eMap->get_acc(p);
}

double target_Info::p_acc(TVector3 p)
{
  return pMap->get_acc(p);
}

double target_Info::pip_acc(TVector3 p)
{
  return pipMap->get_acc(p);
}

double target_Info::getTrans()
{
  return trans;
}

double target_Info::getLum()
{
  return lum;
}

double target_Info::getMass()
{
  return massA;
}

double target_Info::getRadCorr(double Theta, double XB)
{
  if((Theta < 10) || (Theta > 60)){
    std::cout<<"Theta out of bounds.\n";
    return 1;
  }
  if((XB < 0.15) || (XB > 2)){
    std::cout<<"XB out of bounds.\n";
    std::cout<<"XB="<< XB <<"\n";
    return 1;
  }

  double bin_Theta = Theta - 10;
  int b_T = bin_Theta;
  double x_T = bin_Theta - b_T;
  double bin_XB = (XB - 0.15) / 0.05;
  int b_X = bin_XB;
  double x_X = bin_XB - b_X;
  
  double rC = 0;
  rC += (1-x_X) * (1-x_T) * radCorr[b_X][b_T];
  rC += (x_X) * (1-x_T) * radCorr[b_X+1][b_T];
  rC += (1-x_X) * (x_T) * radCorr[b_X][b_T+1];
  rC += (x_X) * (x_T) * radCorr[b_X+1][b_T+1];
  
  return rC;
}


bool target_Info::eVTXInRange(double eVTX)
{
  bool inRange = true;
  if(eVTX < vzMin){
    inRange = false;
  }
  if(eVTX > vzMax){
    inRange = false;
  }
  return inRange;
}

bool target_Info::semiVTXInRange(double eVTX, double leadVTX, double eTheta, double pTheta)
{
  bool inRange = true;
  double diff = eVTX - leadVTX;
 
  if(!eVTXInRange(eVTX)){
    inRange=false;
  }
  if(diff < getVTXMinElectron(eTheta)){
    inRange = false;
  }
  if(diff > getVTXMaxElectron(eTheta)){
    inRange = false;
  }
  if(diff < getVTXMinProton(pTheta)){
    inRange = false;
  }
  if(diff > getVTXMaxProton(pTheta)){
    inRange = false;
  }
  return inRange;
}

bool target_Info::semiVTXInRange(const event_Info myEvent, int leadIndex)
{
  event_Info newEvent = myEvent;
  double eTheta = newEvent.getVector(0).Theta() * 180 / M_PI;
  double pTheta = newEvent.getVector(leadIndex).Theta() * 180 / M_PI;
  return semiVTXInRange(newEvent.getVTX(0),newEvent.getVTX(leadIndex),eTheta,pTheta);
}

bool target_Info::semiFixedVTXInRange(double eVTX, double leadVTX, double eTheta, double pTheta)
{
  bool inRange = true;
  double diff = eVTX - leadVTX;
 
  if(!eVTXInRange(eVTX)){
    inRange=false;
  }
  if(diff < -1){
    inRange = false;
  }
  if(diff > 1){
    inRange = false;
  }
  return inRange;
}

bool target_Info::semiFixedVTXInRange(const event_Info myEvent, int leadIndex)
{
  event_Info newEvent = myEvent;
  double eTheta = newEvent.getVector(0).Theta() * 180 / M_PI;
  double pTheta = newEvent.getVector(leadIndex).Theta() * 180 / M_PI;
  return semiFixedVTXInRange(newEvent.getVTX(0),newEvent.getVTX(leadIndex),eTheta,pTheta);
}

void target_Info::change_vtxMin(double newMin)
{
  vzMin = newMin;
  thick=vzMax-vzMin;
  setLum();
}

void target_Info::change_vtxMax(double newMax)
{
  vzMax = newMax;
  thick=vzMax-vzMin;
  setLum();
}

double target_Info::returnInRangeVTX()
{
  return ((vzMax+vzMin)/2);
}

void target_Info::setLum()
{
  lum = density * thick * ltcc / M_a;
}

void target_Info::fillRadArray()
{
  //Get the file opened
  char radFileName[256]; 
  sprintf(radFileName,"%s/%s_rad_corr.dat",e2adir.c_str(),rad_Name.c_str());
  std::ifstream radFile(radFileName);                                                           
  std::cout<<"Opening Radiation Correctoin File From: \n"<<radFileName<<"\n\n";
  //Check file
  if (!radFile.is_open())
    {
      std::cerr << "Failed to open Radiation Correction file \n"
                << "\n\n Exiting...\n\n";
      exit(-3);
    }

  double Theta,Eprime,Cross,CrossR,Corr,XB;
  //Set to zero and check if it is overwritten
  for(int j = 0; j < 51; j++){
      radCorr[35][j] = 0;
      radCorr[36][j] = 0;
      radCorr[37][j] = 0;
  }
  //Fill with values you find
  while(radFile >> Theta){
    radFile >> Eprime >> Cross >> CrossR >> Corr >> XB;
    if(XB > 2.02){ continue; }
    int binXB = round((XB - 0.15) / 0.05);
    int binTheta = round(Theta - 10);
    radCorr[binXB][binTheta] = Corr;    
  }
  //Some of the files dont have this last bin for some reason
  for(int j = 0; j < 51; j++){
    if(radCorr[35][j] == 0){
      radCorr[35][j] = radCorr[34][j];
    }
    if(radCorr[36][j] == 0){
      radCorr[36][j] = radCorr[35][j];
    }
    if(radCorr[37][j] == 0){
      radCorr[37][j] = radCorr[36][j];
    }
  }
}

double target_Info::getVTXMinElectron(double Theta)
{
  int k = 0;
  while(Theta > eVTX[0][k]){
    k++;
  }
  double x = (Theta-eVTX[0][k])/(eVTX[0][k+1]-eVTX[0][k]);
  double min = ((1-x) * eVTX[1][k]) + (x * eVTX[1][k+1]);
  
  return min;
}

double target_Info::getVTXMaxElectron(double Theta)
{
  int k = 0;
  while(Theta > eVTX[0][k]){
    k++;
  }
  double x = (Theta-eVTX[0][k])/(eVTX[0][k+1]-eVTX[0][k]);
  double max = ((1-x) * eVTX[2][k]) + (x * eVTX[2][k+1]);
  
  return max;
}

double target_Info::getVTXMinProton(double Theta)
{
  int k = 0;
  while(Theta > pVTX[0][k]){
    k++;
  }
  double x = (Theta-pVTX[0][k])/(pVTX[0][k+1]-pVTX[0][k]);
  double min = ((1-x) * pVTX[1][k]) + (x * pVTX[1][k+1]);
  
  return min;
}

double target_Info::getVTXMaxProton(double Theta)
{
  int k = 0;
  while(Theta > pVTX[0][k]){
    k++;
  }
  double x = (Theta-pVTX[0][k])/(pVTX[0][k+1]-pVTX[0][k]);
  double max = ((1-x) * pVTX[2][k]) + (x * pVTX[2][k+1]);
  
  return max;
}

int target_Info::getSecPhi(const double phi){
  
  double nphi=phi;
  if(nphi<=-150) return 0;
  else if(nphi<=-90) return 1;
  else if(nphi<=-30) return 2;
  else if(nphi<=30) return 3;
  else if(nphi<=90) return 4;
  else if(nphi<=150) return 5;
  else return 0; 

}

void target_Info::fillVTXArrayElectron()
{
  //Get the file opened
  char vtxFileName[256]; 
  sprintf(vtxFileName,"%s/%s_et_VTX_cut.txt",e2adir.c_str(),vtx_Name.c_str());
  std::ifstream vtxFile(vtxFileName);                                                           
  std::cout<<"Opening Electron Angle Dependent Vertex Cut File From: \n"<<vtxFileName<<"\n\n";
  //Check file
  if (!vtxFile.is_open())
    {
      std::cerr << "Failed to open Electron Angle Dependent Vertex Cut file \n"
                << "\n\n Exiting...\n\n";
      exit(-3);
    }

  double Theta,vtxMin,vtxMax;
  //Set to zero and check if it is overwritten
  for(int j = 0; j < 50; j++){
    eVTX[0][j]=0;
    eVTX[1][j]=0;
    eVTX[2][j]=0;
  }
  //Fill with values you find
  int k =0;
  while(vtxFile >> Theta){
    vtxFile >> vtxMin >> vtxMax;
    eVTX[0][k]=Theta;
    eVTX[1][k]=vtxMin;
    eVTX[2][k]=vtxMax;
    k++;
  }
  eVTX[0][k]=180;
  eVTX[1][k]=vtxMin;
  eVTX[2][k]=vtxMax;  
}

void target_Info::fillVTXArrayProton()
{
  //Get the file opened
  char vtxFileName[256]; 
  sprintf(vtxFileName,"%s/%s_pt_VTX_cut.txt",e2adir.c_str(),vtx_Name.c_str());
  std::ifstream vtxFile(vtxFileName);                                                           
  std::cout<<"Opening Proton Angle Dependent Vertex Cut File From: \n"<<vtxFileName<<"\n\n";
  //Check file
  if (!vtxFile.is_open())
    {
      std::cerr << "Failed to open Proton Angle Dependent Vertex Cut file \n"
                << "\n\n Exiting...\n\n";
      exit(-3);
    }

  double Theta,vtxMin,vtxMax;
  //Set to zero and check if it is overwritten
  for(int j = 0; j < 50; j++){
    pVTX[0][j]=0;
    pVTX[1][j]=0;
    pVTX[2][j]=0;
  }
  //Fill with values you find
  int k =0;
  while(vtxFile >> Theta){
    vtxFile >> vtxMin >> vtxMax;
    pVTX[0][k]=Theta;
    pVTX[1][k]=vtxMin;
    pVTX[2][k]=vtxMax;
    k++;
  }
  pVTX[0][k]=180;
  pVTX[1][k]=vtxMin;
  pVTX[2][k]=vtxMax;  
}
