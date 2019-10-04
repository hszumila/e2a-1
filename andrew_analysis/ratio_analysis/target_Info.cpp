#include "target_Info.h"
#include "Acceptance.h"
target_Info::target_Info(int A)
{
  
  M_a=(double)A;
  
  if(A==3){
    trans = 0.82;
    vzMax = -0.5;
    vzMin = -2.5;
    density = 0.0655;
    ltcc=10224579.8;
    thick=vzMax-vzMin;
    target = "He3";
  }
  else if(A==4){
    trans = 0.75;
    vzMax = 0.5;
    vzMin = -1.5;
    density = 0.1375;
    ltcc=6895323.02;
    thick=vzMax-vzMin;
    target = "He4";
  }
  else if(A==12){
    trans = 0.53;
    vzMax = 7.0;
    vzMin = 4.0;
    density = 1.786;
    ltcc=17962592.69;
    thick=0.1;
    target = "solid";
  }
  else{
    std::cerr << "The nucleus you have chosen could not be found\n\n Aborting...\n\n";
    exit(-2);
  }
  eMap = new Acceptance(target,4461,2250,"e");
  pMap = new Acceptance(target,4461,2250,"p");
  pipMap = new Acceptance(target,4461,2250,"pip");
  fillRadArray();
  setLum();

}

target_Info::~target_Info()
{
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

bool target_Info::evtxInRange(double eVTX)
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


bool target_Info::vtxInRange(double eVTX, double leadVTX)
{
  bool inRange = true;
  if(eVTX < vzMin){
    inRange = false;
  }
  if(eVTX > vzMax){
    inRange = false;
  }
  if(leadVTX < vzMin){
    inRange = false;
  }
  if(leadVTX > vzMax){
    inRange = false;
  }
  if(abs(eVTX-leadVTX) > 0.5){
    inRange = false;
  }
  return inRange;
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

void target_Info::setLum()
{
  lum = density * thick * ltcc / M_a;
}

void target_Info::fillRadArray()
{
  double Theta,Eprime,Cross,CrossR,Corr,XB;
  radFile.open("");
  for(int i = 0; i < 51; i++){
    for(int j = 0; j < 37; j++){
      radFile >> Theta >> Eprime >> Cross >> CrossR >> Corr >> XB;
      radCorr[j][i] = Corr;
      if(i==0){
	radXB[j] = XB;
      }
      radTheta[i] = Theta;
    }
  }
}
