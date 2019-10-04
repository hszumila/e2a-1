#include "event_Info.h"
#include "Acceptance.h"
srcCut_Info::srcCut_Info()
{

  cutE = false;
  cutLead = false;
  cutRec = false;
  
  doOnlyLeadProtons = false;
  doOnlyLeadNeutrons = false;
  doMinXBCut = false;
  doMaxXBCut = false;
  doMinQSqCut = false;
  doMaxQSqCut = false;
  doMinThetaCut = false;
  doMaxThetaCut = false;
  doMinPoQCut = false;
  doMaxPoQCut = false;
  doMinPMiss = false;
  doMaxPMiss = false;
  doMinMass = false;
  doMaxMass = false;

}

srcCut_Info::~srcCut_Info()
{
}


//Returns false if it does not pass the cut
bool srcCut_Info::applyCut(event_Info &myEvent){
  
}

bool srcCut_Info::passECut(event_Info myEvent)
{
  bool passE = true; 
  if(!passXBCut(myEvent.getXB())){
    passLead=false;
  }
  if(!passQsqCut(myEvent.getQSq())){
    passLead=false;
  }

  return passE;
}

int srcCut_Info::passLeadCut(event_Info myEvent)
{
  int passIndex = -1;
  int numPass = 0;
  for(int i = 1; i < myEvent.getNPar(); i++){
    if(passLeadCutbyIndex(myEvent,i)){
	passIndex = i;
	numPass = 0;
    }
  }
  
  if(numPass < 2){
    return passIndex;
  }
  return -2;

}

bool srcCut_Info::passLeadCutbyIndex(event_Info myEvent, int i)
{
  bool passLead = true;
  if(!passLeadNucleonType(myEvent.getParID(i))){
    passLead=false;
  }
  if(!passThetaCut(myEvent.getThetaPQ(i))){
    passLead=false;
  }
  if(!passPoQCut(myEvent.getPoQ(i))){
    passLead=false;
  }
  if(!passPMissCut(myEvent.getPMiss(i))){
    passLead=false;
  }
  if(!passMassCut(myEvent.getMassMiss(i))){
    passLead=false;
  }
  
  return passE;
}

void srcCut_Info::setOnlyLeadProton()
{
  doOnlyLeadProtons = true;
  cutLead = false;
}

void srcCut_Info::setOnlyLeadNeutron()
{
  doOnlyLeadNeutrons = true;
  cutLead = false;
}

void srcCut_Info::setMinXBCut(double X)
{
  doMinXBCut = true;
  MinXBCut = X;
  cutE = false;
}

void srcCut_Info::setMaxXBCut(double X)
{
  doMaxXBCut = true;
  MaxXBCut = X;
  cutE = false;
}

void srcCut_Info::setMinQSqCut(double X)
{
  doMinQSqCut = true;
  MinQSqCut = X;
  cutE = false;
}

void srcCut_Info::setMaxQSqCut(double X)
{
  doMaxQSqCut = true;
  MaxQSqCut = X;
  cutE = false;
}

void srcCut_Info::setMinThetaCut(double X)
{
  doMinThetaCut = true;
  MinThetaCut = X;
  cutLead = false;
}

void srcCut_Info::setMaxThetaCut(double X)
{
  doMaxThetaCut = true;
  MaxThetaCut = X;
  cutLead = false;
}
void srcCut_Info::setMinPoQCut(double X)
{
  doMinPoQCut = true;
  MinPoQCut = X;
  cutLead = false;
}

void srcCut_Info::setMaxPoQCut(double X)
{
  doMaxPoQCut = true;
  MaxPoQCut = X;
  cutLead = false;
}

void srcCut_Info::setMinPMissCut(double X)
{
  doMinPMissCut = true;
  MinPMissCut = X;
  cutLead = false;
}

void srcCut_Info::setMaxPMissCut(double X)
{
  doMaxPMissCut = true;
  MaxPMissCut = X;
  cutLead = false;
}
void srcCut_Info::setMinMassCut(double X)
{
  doMinMassCut = true;
  MinMassCut = X;
  cutLead = false;
}

void srcCut_Info::setMaxMassCut(double X)
{
  doMaxMassCut = true;
  MaxMassCut = X;
  cutLead = false;
}

bool srcCut_Info::passXBCut(double X)
{
  passCut = true;
  if((doMinXBCut) && (X < MinXBCut)){
    passCut = false;
  }
  if((doMaxXBCut) && (X > MaxXBCut)){
    passCut = false;
  }

  return passCut
}

bool srcCut_Info::passQSqCut(double X)
{
  passCut = true;
  if((doMinQSqCut) && (X < MinQSqCut)){
    passCut = false;
  }
  if((doMaxQSqCut) && (X > MaxQSqCut)){
    passCut = false;
  }

  return passCut
}

bool srcCut_Info::passLeadNucleonType(int ID)
{
  passCut = true;
  if((ID != pcode) && (ID != ncode)){
    passCut = false;
  }
  if(doOnlyLeadProtons && (ID != pcode)){
    passCut = false;
  }
  if(doOnlyLeadNeutrons && (ID != ncode)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passThetaCut(double X)
{
  passCut = true;
  if((doMinThetaCut) && (X < MinThetaCut)){
    passCut = false;
  }
  if((doMaxThetaCut) && (X > MaxThetaCut)){
    passCut = false;
  }

  return passCut
}

bool srcCut_Info::passPoQCut(double X)
{
  passCut = true;
  if((doMinPoQCut) && (X < MinPoQCut)){
    passCut = false;
  }
  if((doMaxPoQCut) && (X > MaxPoQCut)){
    passCut = false;
  }

  return passCut
}

bool srcCut_Info::passPMissCut(double X)
{
  passCut = true;
  if((doMinPMissCut) && (X < MinPMissCut)){
    passCut = false;
  }
  if((doMaxPMissCut) && (X > MaxPMissCut)){
    passCut = false;
  }

  return passCut
}

bool srcCut_Info::passMassCut(double X)
{
  passCut = true;
  if((doMinMassCut) && (X < MinMassCut)){
    passCut = false;
  }
  if((doMaxMassCut) && (X > MaxMassCut)){
    passCut = false;
  }

  return passCut
}


