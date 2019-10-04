#ifndef __PARTINFO_H__
#define __PARTINFO_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdlib>

#include "TVector3.h"

#include "e2a_constants.h"

class part_Info
{
 public:
  part_Info(int particleID, double mom_x, double mom_y, double mom_z, double vertex);
  ~part_Info();
  //Functions to test the type of particle
  bool isElectron();
  bool isNucleon();
  bool isProton();
  bool isNeutron();
  bool isPion();
  bool isDelta();
  //Functions to get values
  int getParID();
  double getPX();
  double getPY();
  double getPZ();
  double getVTX();
  TVector3 getVector();

 
 private:
  int parID;
  TVector3 mom;
  double vtx;

};

#endif
