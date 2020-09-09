#include "srcCut_Info.h"
#include "e2a_constants.h"
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
  doMinPMissCut = false;
  doMaxPMissCut = false;
  doMinMassCut = false;
  doMaxMassCut = false;

}

srcCut_Info::~srcCut_Info()
{
}


//Returns false if it does not pass the cut
int srcCut_Info::passCutReorder(event_Info &myEvent)
{  
  //If it is only a cut on the electrons
  if(cutE && !cutLead){
    if(passECut(myEvent)){
      myEvent.clearNonElectron();
      return 1;
    }
    return 0;
  }    //If it is a cut on the electron and the lead
  else if(cutE && cutLead){
    if(!passECut(myEvent)){
	return 0;
      }
    else{
      int leadIndex = passLeadCut(myEvent);      
      if(leadIndex > 0){
	myEvent.setLeadandClear(leadIndex);	
	return 1;
      }
      else if(leadIndex == -2){
	return 2;
      }
      else{
	return 0;
      }      
    }
  }
}

int srcCut_Info::passCutConst(const event_Info myEvent)
{
  //If it is only a cut on the electrons
  if(cutE && !cutLead){
    if(passECut(myEvent)){
      return 1;
    }
    return 0;
  }    //If it is a cut on the electron and the lead
  else if(cutE && cutLead){
   
    if(!passECut(myEvent)){
	return 0;
      }
    else{
      int leadIndex = passLeadCut(myEvent);
      
      if(leadIndex > 0){
	return 1;
      }
      else if(leadIndex == -2){
	return 2;
      }
      else{
	return 0;
      }      
    }
  }

}


void srcCut_Info::makeInclCut()
{
  setMinQSqCut(1.4);
  setMaxQSqCut(5);
}

void srcCut_Info::makeStandardSemiCut()
{
  setMinXBCut(1.2);
  setMaxXBCut(2);
  setMaxThetaCut(25);
  setMinPoQCut(0.62);
  setMaxPoQCut(0.96);
  setMinPMissCut(0.3);
  setMaxPMissCut(0.6);
  setMinMassCut(0.9);
  setMaxMassCut(1.1);
}

void srcCut_Info::makeNewSemiCut()
{
  setMaxXBCut(2);
  setMaxThetaCut(25);
  setMinPoQCut(0.62);
  setMinPMissCut(0.3);
  setMaxPMissCut(0.6);
  setMinMassCut(0.9);
  setMaxMassCut(1.1);
}

void srcCut_Info::makeNewSemiCutNoPMiss()
{
  setMaxXBCut(2);
  setMaxThetaCut(25);
  setMinPoQCut(0.62);
  setMaxPMissCut(2);
  setMinMassCut(0.9);
  setMaxMassCut(1.1);
}

void srcCut_Info::makeGausCut()
{
  setMinXBCut(0.15);
  setMaxXBCut(2);

  TRandom3 myRand(0);
  /*  setMaxThetaCut(myRand.Gaus(25,5));
  setMinPoQCut(myRand.Gaus(0.62,0.05));
  setMinPMissCut(myRand.Gaus(0.3,0.025));
  setMaxPMissCut(myRand.Gaus(0.6,0.025));
  setMinMassCut(myRand.Gaus(0.9,0.05));
  setMaxMassCut(myRand.Gaus(1.1,0.05));
  */

  setMaxThetaCut(myRand.Uniform(20,30));
  setMinThetaCut(myRand.Uniform(0,5));
  setMinPoQCut(myRand.Uniform(0.57,0.67));
  setMinPMissCut(myRand.Uniform(0.275,0.325));
  setMaxPMissCut(myRand.Uniform(0.575,0.625));
  setMinMassCut(myRand.Uniform(0.85,0.95));
  setMaxMassCut(myRand.Uniform(1.05,1.15));
}

void srcCut_Info::makeLightCut()
{
  setMaxXBCut(2);
  setMaxThetaCut(25);
  setMinPoQCut(0.62);
  setMaxPMissCut(2);
}

void srcCut_Info::makeEG2Cut()
{
  setMinPMissCut(0.4);
  setMaxPMissCut(0.6);
  setMinXBCut(1.2);
  setMaxThetaCut(25);
  setMinPoQCut(0.62);
  setMaxPoQCut(0.96);
  setMaxMassCut(1.1);
}

bool srcCut_Info::passECut(event_Info myEvent)
{
  bool passE = true; 
  if(!passXBCut(myEvent.getXB())){
    passE=false;
  }
  if(!passQSqCut(myEvent.getQSq())){
    passE=false;
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
  
  return passLead;
}

void srcCut_Info::setOnlyLeadProton()
{
  doOnlyLeadProtons = true;
  cutLead = true;
}

void srcCut_Info::setOnlyLeadNeutron()
{
  doOnlyLeadNeutrons = true;
  cutLead = true;
}

void srcCut_Info::setMinXBCut(double X)
{
  doMinXBCut = true;
  MinXBCut = X;
  cutE = true;
}

void srcCut_Info::setMaxXBCut(double X)
{
  doMaxXBCut = true;
  MaxXBCut = X;
  cutE = true;
}

void srcCut_Info::setMinQSqCut(double X)
{
  doMinQSqCut = true;
  MinQSqCut = X;
  cutE = true;
}

void srcCut_Info::setMaxQSqCut(double X)
{
  doMaxQSqCut = true;
  MaxQSqCut = X;
  cutE = true;
}

void srcCut_Info::setMinThetaCut(double X)
{
  doMinThetaCut = true;
  MinThetaCut = X;
  cutLead = true;
}

void srcCut_Info::setMaxThetaCut(double X)
{
  doMaxThetaCut = true;
  MaxThetaCut = X;
  cutLead = true;
}
void srcCut_Info::setMinPoQCut(double X)
{
  doMinPoQCut = true;
  MinPoQCut = X;
  cutLead = true;
}

void srcCut_Info::setMaxPoQCut(double X)
{
  doMaxPoQCut = true;
  MaxPoQCut = X;
  cutLead = true;
}

void srcCut_Info::setMinPMissCut(double X)
{
  doMinPMissCut = true;
  MinPMissCut = X;
  cutLead = true;
}

void srcCut_Info::setMaxPMissCut(double X)
{
  doMaxPMissCut = true;
  MaxPMissCut = X;
  cutLead = true;
}
void srcCut_Info::setMinMassCut(double X)
{
  doMinMassCut = true;
  MinMassCut = X;
  cutLead = true;
}

void srcCut_Info::setMaxMassCut(double X)
{
  doMaxMassCut = true;
  MaxMassCut = X;
  cutLead = true;
}

//Functions to get the cut values
double srcCut_Info::getMinXB()
{
  if(doMinXBCut == true){return MinXBCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxXB()
{
  if(doMaxXBCut == true){return MaxXBCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMinQSq()
{
  if(doMinQSqCut == true){return MinQSqCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxQSq()
{
  if(doMaxQSqCut == true){return MaxQSqCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}
double srcCut_Info::getMinTheta()
{
  if(doMinThetaCut == true){return MinThetaCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxTheta()
{
  if(doMaxThetaCut == true){return MaxThetaCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}
double srcCut_Info::getMinPoQ()
{
  if(doMinPoQCut == true){return MinPoQCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxPoQ()
{
  if(doMaxPoQCut == true){return MaxPoQCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}
double srcCut_Info::getMinPMiss()
{
  if(doMinPMissCut == true){return MinPMissCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxPMiss()
{
  if(doMaxPMissCut == true){return MaxPMissCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}
double srcCut_Info::getMinMass()
{
  if(doMinMassCut == true){return MinMassCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}

double srcCut_Info::getMaxMass()
{
  if(doMaxMassCut == true){return MaxMassCut;}
  std::cerr<<"Attempting to retrieve a cut value that was not set.\n Returning 0.\n\n";
  return 0;
}
//Check to see if the event passes the cut
bool srcCut_Info::passXBCut(double X)
{
  bool passCut = true;
  if((doMinXBCut) && (X < MinXBCut)){
    passCut = false;
  }
  if((doMaxXBCut) && (X > MaxXBCut)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passQSqCut(double X)
{
  bool passCut = true;
  if((doMinQSqCut) && (X < MinQSqCut)){
    passCut = false;
  }
  if((doMaxQSqCut) && (X > MaxQSqCut)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passLeadNucleonType(int ID)
{
  bool passCut = true;
  if((ID != pCode) && (ID != nCode)){
    passCut = false;
  }
  if(doOnlyLeadProtons && (ID != pCode)){
    passCut = false;
  }
  if(doOnlyLeadNeutrons && (ID != nCode)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passThetaCut(double X)
{
  bool passCut = true;
  if((doMinThetaCut) && (X < MinThetaCut)){
    passCut = false;
  }
  if((doMaxThetaCut) && (X > MaxThetaCut)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passPoQCut(double X)
{
  bool passCut = true;
  if((doMinPoQCut) && (X < MinPoQCut)){
    passCut = false;
  }
  if((doMaxPoQCut) && (X > MaxPoQCut)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passPMissCut(double X)
{
  bool passCut = true;
  if((doMinPMissCut) && (X < MinPMissCut)){
    passCut = false;
  }
  if((doMaxPMissCut) && (X > MaxPMissCut)){
    passCut = false;
  }

  return passCut;
}

bool srcCut_Info::passMassCut(double X)
{
  bool passCut = true;
  if((doMinMassCut) && (X < MinMassCut)){
    passCut = false;
  }
  if((doMaxMassCut) && (X > MaxMassCut)){
    passCut = false;
  }

  return passCut;
}


