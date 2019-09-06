#ifndef __EVENTINFO_H__
#define __EVENTINFO_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdlib>

#include "TVector3.h"

#include "e2a_constants.h"

class event_Info
{
 public:
  event_Info(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19], double minimumXb, double minimumPMiss, bool onlyAcceptLeadProtons, bool onlyAcceptLeadNeutrons);
  //event_Info(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19]);
  ~event_Info();
  void fillValues(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19], double minimumXb, double minimumPMiss, bool onlyAcceptLeadProtons, bool onlyAcceptLeadNeutrons);
  //Functions to test the type of particle
  bool isNucleon(int j);
  bool isProton(int j);
  bool isNeutron(int j);
  bool isPion(int j);
  bool isDelta(int j);
  int getDeltaType(int j, int k);
  //Functions to get values
  int getNumDeltas();
  int getNPar();
  int getParID(int i);
  double getPX(int i);
  double getPY(int i);
  double getPZ(int i);
  double getVTX(int i);
  //Functions to related to lead and recoil
  void setLead(int i);
  void setRec(int i);
  int getWhichLead();
  bool isLeadbyIndex(int i);
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
  int parID[19];
  double xB;
  double px[19];
  double py[19];
  double pz[19];
  double vtx[19];

  double minP;
  double minX;
  bool onlyLeadProtons;
  bool onlyLeadNeutrons;

  bool vtxMatch(int j, int k);
  bool checkDeltaWithMass(const double dMass);
  //Functions for manipulating the set
  void combineParticle(int j, int k, int newParID);
  void moveEntryForward(int startIndex, int endIndex);
  void copy(int indexOverWrite, int indexCopy);
  void changeNPar(int new_nPar);
  void mergParInArrays(int j, int k, int newParID, TVector3 vDelta, double new_vtx);
};

#endif
