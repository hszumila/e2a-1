#ifndef __TARGETINFO_H__
#define __TARGETINFO_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdlib>

#include "TFile.h"
#include "TH3.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "Acceptance.h"
#include "e2a_constants.h"

class TH3D;
class TFile;

const int nBinsMom=100;
const int nBinsCos=200;
const int nBinsPhi=360;

class target_Info
{
 public:
  target_Info(int A);
  ~target_Info();
  double incl_acc(const TVector3 ve);
  double semi_acc(const TVector3 ve, const TVector3 vLead);
  double e_acc(TVector3 p);
  double p_acc(TVector3 p);
  double pip_acc(TVector3 p);
  double getTrans();
  double getLum();
  double getRadCorr();
  bool evtxInRange(double eVTX);
  bool vtxInRange(double eVTX, double leadVTX);
  void change_vtxMin(double newMin);
  void change_vtxMax(double newMax);
  
  
 private:
  double trans;
  double lum;
  double density;
  double M_a;
  double ltcc; //live-time-corrected-charge
  double vzMax;
  double vzMin;
  double thick;
  std::string target;
  std::ifstream radFile;
  double radXB[37];
  double radTheta[51];
  double radCorr[37][51];
  Acceptance * eMap = NULL; 
  Acceptance * pMap = NULL;
  Acceptance * pipMap = NULL;
  void setLum();
  void fillRadArray();
  
};

#endif
