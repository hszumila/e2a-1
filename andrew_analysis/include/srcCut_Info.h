#ifndef __SRCCUTINFO_H__
#define __SRCCUTINFO_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdlib>

#include "TVector3.h"
#include "TRandom3.h"

#include "e2a_constants.h"
#include "event_Info.h"

class srcCut_Info
{
 public:
  srcCut_Info();
  ~srcCut_Info();
  int passCutReorder(event_Info &myEvent);
  int passCutConst(const event_Info myEvent);
  //Predefined Cuts
  void makeInclCut();
  void makeStandardSemiCut();
  void makeNewSemiCut();
  void makeGausCut();
  void makeLightCut();
  //Check what it passed
  bool passECut(event_Info myEvent);
  int passLeadCut(event_Info myEvent);
  bool passLeadCutbyIndex(event_Info myEvent, int i);
  //Functions to set a cut
  void setOnlyLeadProton();
  void setOnlyLeadNeutron();
  void setMinXBCut(double X);
  void setMaxXBCut(double X);
  void setMinQSqCut(double X);
  void setMaxQSqCut(double X);
  void setMinThetaCut(double X);
  void setMaxThetaCut(double X);
  void setMinPoQCut(double X);
  void setMaxPoQCut(double X);
  void setMinPMissCut(double X);
  void setMaxPMissCut(double X);
  void setMinMassCut(double X);
  void setMaxMassCut(double X);
 
 private:
  bool cutE;
  bool cutLead;
  bool cutRec;
  
  bool doOnlyLeadProtons;
  bool doOnlyLeadNeutrons;
  bool doMinXBCut;
  bool doMaxXBCut;
  bool doMinQSqCut;
  bool doMaxQSqCut;
  bool doMinThetaCut;
  bool doMaxThetaCut;
  bool doMinPoQCut;
  bool doMaxPoQCut;
  bool doMinPMissCut;
  bool doMaxPMissCut;
  bool doMinMassCut;
  bool doMaxMassCut;

  double MinXBCut;
  double MaxXBCut;
  double MinQSqCut;
  double MaxQSqCut;
  double MinThetaCut;
  double MaxThetaCut;
  double MinPoQCut;
  double MaxPoQCut;
  double MinPMissCut;
  double MaxPMissCut;
  double MinMassCut;
  double MaxMassCut;

  bool passLeadNucleonType(int ID);
  bool passXBCut(double X);
  bool passQSqCut(double X);
  bool passThetaCut(double X);
  bool passPoQCut(double X);
  bool passPMissCut(double X);
  bool passMassCut(double X);
  
};

#endif
