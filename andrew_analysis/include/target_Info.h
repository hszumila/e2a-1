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
#include "Fiducial.h"
#include "e2a_constants.h"
#include "event_Info.h"

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
  bool pass_incl_fid(const TVector3 ve);
  bool pass_semi_fid(const TVector3 ve, const TVector3 vLead);
  double e_acc(TVector3 p);
  double p_acc(TVector3 p);
  double pip_acc(TVector3 p);
  double getTrans();
  double getLum();
  double getMass();
  double getRadCorr(double Theta, double XB);

  bool eVTXInRange(double eVTX);
  bool semiVTXInRange(double eVTX, double leadVTX, double eTheta, double pTheta);
  bool semiVTXInRange(const event_Info myEvent, int leadIndex);
  bool semiFixedVTXInRange(double eVTX, double leadVTX, double eTheta, double pTheta);
  bool semiFixedVTXInRange(const event_Info myEvent, int leadIndex);
  void change_vtxMin(double newMin);
  void change_vtxMax(double newMax);
  double returnInRangeVTX();
  
  int getSecPhi(const double phi);
  
 private:
  double trans;
  double lum;
  double density;
  double M_a;
  double massA;
  double ltcc; //live-time-corrected-charge
  double vzMax;
  double vzMin;
  double thick;
  std::string e2adir;
  std::string acc_Name;
  std::string fid_Name;
  std::string rad_Name;
  std::string vtx_Name;
  double radCorr[38][51];
  double eVTX[3][50];
  double pVTX[3][50];
  Acceptance * eMap = NULL; 
  Acceptance * pMap = NULL;
  Acceptance * pipMap = NULL;
  Fiducial * targFid = NULL;
  void setLum();
  std::ifstream radFile;
  void fillRadArray();

  double getVTXMinElectron(double Theta);
  double getVTXMaxElectron(double Theta);
  double getVTXMinProton(double Theta);
  double getVTXMaxProton(double Theta);
  void fillVTXArrayElectron();
  void fillVTXArrayProton();
  
};

#endif
