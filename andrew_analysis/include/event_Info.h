#ifndef __EVENTINFO_H__
#define __EVENTINFO_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "TVector3.h"

#include "e2a_constants.h"
#include "part_Info.h"

class event_Info
{
 public:
  event_Info(int numberOfParticles, int particleIDs[19], double Xb, double Q2, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19]);
  ~event_Info();
  //Functions to test the type of particle
  bool isElectron(int i);
  bool isNucleon(int i);
  bool isProton(int i);
  bool isNeutron(int i);
  bool isPion(int i);
  bool isDelta(int i);
  int getDeltaType(int j, int k);
  //Functions to get values
  int getNumDeltas();
  int getNPar();
  int getParID(int i);
  double getPX(int i);
  double getPY(int i);
  double getPZ(int i);
  TVector3 getVector(int i);
  double getVTX(int i);
  double getThetaPQ(int i);
  double getPoQ(int i);
  double getPMiss(int i);
  double getMassMiss(int i);
  double getXB();
  double getQSq();
  //Functions to related to lead and recoil
  void clearNonElectron();
  void setLead(int i);
  void setLeadandClear(int i);
  void setRec(int i);
  //Functions related to the delta
  void findAndMergeDeltas();
  void addDelta(int j, int k);
  bool doesNucleonMatchPion(int i);
  int getWhichPionMatch(int i);
  double getMassNucleonPion(int j, int k);
  bool checkDeltaWithIndex(int j, int k);

 
 private:
  int nPar;
  int nDeltas;
  double xB;
  double QSq;
  std::vector<part_Info> parList;

  bool vtxMatch(int j, int k);
  bool checkDeltaWithMass(const double dMass);
  //Functions for manipulating the set
  void combineParticle(int j, int k, int newParID);
  void moveEntry(int startIndex, int endIndex);
  void clearAbove(int startIndex);
  void copy(int indexOverWrite, int indexCopy);
  void changeNPar(int new_nPar);
  void mergeParInArrays(int j, int k, part_Info merged);
};

#endif
