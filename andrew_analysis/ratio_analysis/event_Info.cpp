#include "event_Info.h"
#include "Acceptance.h"
event_Info::event_Info(int numberOfParticles, int particleIDs[19], double Xb, double Q2, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19])
{
  nPar=numberOfParticles;
  nDeltas = 0;
  xB=Xb;
  QSq=Q2;
  for(int i = 0; i < nPar; i++){
    part_Info thisPart(particleIDs[i], mom_x[i], mom_y[i], mom_z[i], vertecies[i]);
    parList.push_back(thisPart);
  }

}

event_Info::~event_Info()
{
}

//Particle ID functions
bool event_Info::isElectron(int i)
{
  return parList.at(i).isElectron();
}

bool event_Info::isNucleon(int i)
{
  return parList.at(i).isNucleon();
}

bool event_Info::isProton(int i)
{
  return parList.at(i).isProton();  
}

bool event_Info::isNeutron(int i)
{
  return parList.at(i).isNeutron();
}


bool event_Info::isPion(int i)
{
return parList.at(i).isPion();
}

bool event_Info::isDelta(int i)
{
return parList.at(i).isDelta();
}

int event_Info::getDeltaType(int j, int k){
  int nucleonID = parList.at(j).getParID();
  int pionID = parList.at(j).getParID();
  
  if((nucleonID==pCode) && (pionID==pipCode) ){
    return dppCode;
  }
  if((nucleonID==pCode) && (pionID==pi0Code) ){
    return dpCode;
  }
  if((nucleonID==pCode) && (pionID==pimCode) ){
    return d0Code;
  }
  if((nucleonID==nCode) && (pionID==pipCode) ){
    return dpCode;
  }
  if((nucleonID==nCode) && (pionID==pi0Code) ){
    return d0Code;
  }
  if((nucleonID==nCode) && (pionID==pimCode) ){
    return dmCode;
  }

  std::cerr<<"Got a delta with incorrect codes, returning 0.\n";

  return 0;
  
}


//Functions to get numbers
int event_Info::getNumDeltas()
{
  return nDeltas;
}

int event_Info::getNPar()
{
  return nPar;
}

int event_Info::getParID(int i)
{
  return parList.at(i).getParID();
}

double event_Info::getPX(int i)
{
  return parList.at(i).getPX();
}

double event_Info::getPY(int i)
{
  return parList.at(i).getPY();
}

double event_Info::getPZ(int i)
{
  return parList.at(i).getPZ();
}

TVector3 event_Info::getVector(int i)
{
  return parList.at(i).getVector();
}

double event_Info::getVTX(int i)
{
  return parList.at(i).getVTX();
}

double event_Info::getThetaPQ(int i)
{
  TVector3 vLead = getVector(i);
  TVector3 vq = vBeam - getVector(0);
  double ThetaPQ = vq.Angle(vLead) * (180/M_PI);
  
  return ThetaPQ;
}

double event_Info::getPoQ(int i)
{
  TVector3 vLead = getVector(i);
  TVector3 vq = vBeam - getVector(0);
  double PoQ = vLead.Mag()/vq.Mag();
  
  return PoQ;
}

double event_Info::getPMiss(int i)
{
  TVector3 vLead = getVector(i);
  TVector3 vq = vBeam - getVector(0);
  TVector3 vMiss = vLead - vq;
  double PMiss = vMiss.Mag();

  return PMiss;
}

double event_Info::getMassMiss(int i)
{ 
  TVector3 vLead = getVector(i);
  TVector3 vq = vBeam - getVector(0);
  double omega = vBeam.Mag() - getVector(0).Mag();
  TVector3 vMiss = vLead - vq;
  double eLead = sqrt(vLead.Mag2() + (mP*mP));
  double eMiss = omega + mP + mN - eLead;
  double mMiss = sqrt((eMiss*eMiss)-vMiss.Mag2());

  return mMiss;
}

double event_Info::getXB()
{
  return xB;
}

double event_Info::getQSq()
{
  return QSq;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions to find and set the lead
void event_Info::clearNonElectron()
{
  clearAbove(1);
}

void event_Info::setLead(int i)
{
  moveEntry(i,1);
}

void event_Info::setLeadandClear(int i)
{
  moveEntry(i,1);
  clearAbove(2);
}

void event_Info::setRec(int i)
{
  moveEntry(i,2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions to add the delta
//This functions takes the nucleons in the order that they came and looks
//For a pion that matches to a delta with it. If it is found, then they are
//Matched and the search continues. If two pions match a nucleon, then the
//Match is not made. If two nucleons match with a pion, then the first pion
//To match is joined with the pion
void event_Info::findAndMergeDeltas()
{
  for(int i = 1; i < nPar; i++){
    int l = getWhichPionMatch(i);
    if(l > 0){
      addDelta(i,l);
    }
  }
}

//Adds a delta to the end of the list of particles by merging entry j and k
void event_Info::addDelta(int j, int k)
{
  int dType = getDeltaType(j,k);
  combineParticle(j,k,dType);
  nDeltas++;
}
//Returns true if the nucleon can be matched with a pion to create a delta
bool event_Info::doesNucleonMatchPion(int i){
  if(getWhichPionMatch(i) > 0 ){
    return true;
  }
  return false;
}
//You input the index of a test nucleon. This function looks through all indecies to see
//if it matches with a pion. It returns the index if it matches with a pion. It returns
//-1 if there is no match. It returns -2 if more than one match is made.
int event_Info::getWhichPionMatch(int i){
  int numMatches = 0;
  int MatchIndex = -1;
  for(int l = 1; l < nPar; l++){
    if(checkDeltaWithIndex(i,l)){
      MatchIndex = l;
      numMatches++;
    }
  }
  if(numMatches < 2){
    return MatchIndex;
  }
  return -2;
}
//Returns the invariant mass of a nucleon at index j and a pion at index k
double event_Info::getMassNucleonPion(int j, int k){

  TVector3 vN = parList.at(j).getVector();
  TVector3 vpi = parList.at(k).getVector();
  TVector3 vD_test = vN + vpi;
  double eN = sqrt(vN.Mag2()+(mN*mN));
  double epi = sqrt(vpi.Mag2()+(mpc*mpc));
  double eD_test = eN + epi;
  double mD_test = sqrt((eD_test*eD_test)-vD_test.Mag2());
  
  return mD_test;    
}
//Checks to see if a nucleon(j) and pion(k) could be considered a delta
bool event_Info::checkDeltaWithIndex(int j, int k){
  if(vtxMatch(j,k)){
    if(isNucleon(j) && isPion(k)){
      return checkDeltaWithMass(getMassNucleonPion(j,k));
    }
  }
  return false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Private Functions
bool event_Info::vtxMatch(int j, int k)
{
  if(abs(parList.at(j).getVTX()-parList.at(k).getVTX()) < 0.5){
    return true;
  }

  return false;
  
}

bool event_Info::checkDeltaWithMass(const double dMass){
if((dMass > (mD-wD)) && (dMass < (mD+wD))){
      return true;
    }
  return false;
}

void event_Info::combineParticle(int j, int k, int newParID)
{
  changeNPar(nPar-1);
  TVector3 vDelta = parList.at(j).getVector() + parList.at(k).getVector();
  double new_vtx = (parList.at(j).getVTX()+parList.at(k).getVTX())/2;
  part_Info merged(newParID, vDelta.X(), vDelta.X(), vDelta.X(), new_vtx);
  mergeParInArrays(j,k,merged);
}

void event_Info::moveEntry(int startIndex, int endIndex)
{
  part_Info temp = parList.at(startIndex);
  parList.erase(parList.begin()+startIndex);
  parList.insert(parList.begin()+endIndex,temp);

}

void event_Info::clearAbove(int startIndex)
{
  for(int i = startIndex; i < nPar; i++){
    parList.erase(parList.begin()+startIndex);
  }
  nPar = startIndex;
}

void event_Info::changeNPar(int new_nPar)
{
  nPar = new_nPar;
}

void event_Info::mergeParInArrays(int j, int k, part_Info merged)
{
  parList.erase(parList.begin()+j);
  parList.erase(parList.begin()+k);
  parList.push_back(merged);
}
