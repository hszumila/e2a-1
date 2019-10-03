#include "event_Info.h"
#include "Acceptance.h"
srcCut_Info::srcCut_Info()
{
  doOnlyLeadProtons = false;
  doOnlyLeadNeutrons = false;
  doMinXBCut = false;
  doMaxXBCut = false;
  doMaxThetaCut = false;
  doMinPoQCut = false;
  doMaxPoQCut = false;
  do
}

srcCut_Info::~srcCut_Info()
{
}


//Finds out which particle is a lead nucleon in the entire event
//If none are lead return -1. If more than one is lead return -2.
int event_Info::getWhichLead(){
  int numLead = 0;
  int leadIndex = -1;
  for(int j = 1; j < nPar; j++){
    if(isLeadbyIndex(j)){
      leadIndex = j;
      numLead++;
    }
  }
  if(numLead < 2){
    return leadIndex;
  }
  return -2;
}
//Tests to see if the nucleon at index i is a lead
bool event_Info::isLeadbyIndex(int i){
  
  const TVector3 vBeam(0.,0.,4.461);
  TVector3 ve(px[0],py[0],pz[0]);
  TVector3 vLead(px[i],py[i],pz[i]);
  TVector3 vq = vBeam - ve;
  TVector3 vMiss = vLead - vq;
  bool passLead = true;
  
  if(!isNucleon(i)){
    passLead=false;
  }
  if(onlyLeadProtons && (!isProton(i))){
    passLead=false;
  }
  if(onlyLeadNeutrons && (!isNeutron(i))){
    passLead=false;
  }
  if(xB < minX){
    passLead=false;
  }
  if(xB > 2){
  passLead=false;
  }
  if(vq.Angle(vLead)>((M_PI*25)/180)){
    passLead=false;
  }
  if((vLead.Mag()/vq.Mag())<(0.62)){
    passLead=false;
  }
  //if((vLead.Mag()/vq.Mag())>(0.96)){
  //  passLead=false;
  //}
  if(vMiss.Mag() < (minP)){
    passLead=false;
  }
  //if(passLead){std::cout<<passLead<<"\n";}
  if(vMiss.Mag() > (0.6)){
    passLead=false;
  }
  return passLead;
}

