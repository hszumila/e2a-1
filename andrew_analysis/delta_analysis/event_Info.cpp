#include "event_Info.h"
#include "Acceptance.h"
event_Info::event_Info(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19], double minimumXb, double minimumPMiss, bool onlyAcceptLeadProtons, bool onlyAcceptLeadNeutrons)
{
  fillValues(numberOfParticles, particleIDs, Xb, mom_x, mom_y, mom_z, vertecies, minimumXb, minimumPMiss, onlyAcceptLeadProtons, onlyAcceptLeadNeutrons);
}

//event_Info::event_Info(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19])
//{
//fillValues(numberOfParticles, particleIDs, Xb, mom_x, mom_y, mom_z, vertecies, 1.2, 0.3, true, false);
//}

event_Info::~event_Info()
{
}

void event_Info::fillValues(int numberOfParticles, int particleIDs[19], double Xb, double mom_x[19], double mom_y[19], double mom_z[19], double vertecies[19], double minimumXb, double minimumPMiss, bool onlyAcceptLeadProtons, bool onlyAcceptLeadNeutrons)
{
  nPar=numberOfParticles;
  nDeltas = 0;
  xB=Xb;
  for(int i = 0; i < 19; i++){
    parID[i]=particleIDs[i];
    px[i]=mom_x[i];
    py[i]=mom_y[i];
    pz[i]=mom_z[i];
    vtx[i]=vertecies[i];
  }

  minX = minimumXb;
  minP = minimumPMiss;
  onlyLeadProtons = onlyAcceptLeadProtons;
  onlyLeadNeutrons = onlyAcceptLeadNeutrons;

}

//Particle ID functions
bool event_Info::isNucleon(int j){
  int ID = parID[j];
  if(ID == pCode){
    return true;
  }
  else if(ID == nCode){
    return true;
  }
  return false;
}

bool event_Info::isProton(int j){
  int ID = parID[j];
  if(ID == pCode){
    return true;
  }
  return false;
}

bool event_Info::isNeutron(int j){
  int ID = parID[j];
  if(ID == nCode){
    return true;
  }
  return false;
}


bool event_Info::isPion(int j){
  int ID = parID[j];
  if(ID == pipCode){
    return true;
  }
  else if(ID == pimCode){
    return true;
  }
  else if(ID == pi0Code){
    return true;
  }
  return false;
}

bool event_Info::isDelta(int j){
  int ID = parID[j];
  if(ID == dppCode){
    return true;
  }
  else if(ID == dpCode){
    return true;
  }
  else if(ID == d0Code){
    return true;
  }
  else if(ID == dmCode){
    return true;
  }
  return false;
}

int event_Info::getDeltaType(int j, int k){
  int nucleonID = parID[j];
  int pionID = parID[k];
  
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
  return parID[i];
}

double event_Info::getPX(int i)
{
  return px[i];
}

double event_Info::getPY(int i)
{
  return py[i];
}

double event_Info::getPZ(int i)
{
  return pz[i];
}

double event_Info::getVTX(int i)
{
  return vtx[i];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions to find and set the lead
void event_Info::setLead(int i)
{
  moveEntryForward(i,1);
}

void event_Info::setRec(int i)
{
  moveEntryForward(i,2);
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

  TVector3 vN(px[j],py[j],pz[j]);
  TVector3 vpi(px[j],py[j],pz[j]);
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
  if(abs(vtx[j]-vtx[k]) < 0.5){
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
  TVector3 vDelta((px[j]+px[k]),(py[j]+py[k]),(pz[j]+pz[k]));
  double new_vtx = (vtx[j]+vtx[k])/2;
  mergParInArrays(j,k,newParID,vDelta,new_vtx);
}

void event_Info::moveEntryForward(int startIndex, int endIndex)
{
  if(endIndex > startIndex){
    std::cerr<<"You are trying to move a particle forward in the index from index "<<startIndex<<" to index "<<endIndex<<".\n Aborting..";
    exit(-2);
  }

  int startID = parID[startIndex];
  double startpx = px[startIndex];
  double startpy = py[startIndex];
  double startpz = pz[startIndex];
  double startvtx = vtx[startIndex];

  for(int l = startIndex; l > endIndex; l--){
    copy(l,(l-1));
  }

  parID[endIndex]=startID;
  px[endIndex]=startpx;
  py[endIndex]=startpy;
  pz[endIndex]=startpz;;
  vtx[endIndex]=startvtx;  

}

void event_Info::copy(int indexOverWrite, int indexCopy){
    parID[indexOverWrite]=parID[indexCopy];
    px[indexOverWrite]=px[indexCopy];
    py[indexOverWrite]=py[indexCopy];
    pz[indexOverWrite]=pz[indexCopy];
    vtx[indexOverWrite]=vtx[indexCopy];  
}
  
void event_Info::changeNPar(int new_nPar)
{
  nPar = new_nPar;
}

void event_Info::mergParInArrays(int j, int k, int newParID, TVector3 vDelta, double new_vtx)
{
  //Check which particle comes first in the list
  int a,b;
  if(j<k){
    a = j;
    b = k;
  }
  else if(k<j){
    a = k;
    b = j;
  }
  else{
    std::cerr<<"The program has picked the same particle to combine with itself.\n Aborting..";
    exit(-2);
  } 

  //Now I move all particles between the two up by one
  for(int i = a; i < (b-1) ; i++){
      parID[i]=parID[i+1];
      px[i]=px[i+1];
      py[i]=py[i+1];
      pz[i]=pz[i+1];      
      vtx[i]=vtx[i+1];
    }
  //Now I move all particles after the two up by two
  for(int i = (b-1); i < (nPar-1) ; i++){
      parID[i]=parID[i+2];
      px[i]=px[i+2];
      py[i]=py[i+2];
      pz[i]=pz[i+2];
      vtx[i]=vtx[i+2];
  }
  //Now set anything else equal to zero
  for(int i = (nPar-1); i < 19; i++){
    parID[i]=0;
    px[i]=0;
    py[i]=0;
    pz[i]=0;
    vtx[i]=0;

  }
  //Finally, make the last particle the delta
  parID[nPar-1] = newParID;
  px[nPar-1]=vDelta.X();
  py[nPar-1]=vDelta.Y();
  pz[nPar-1]=vDelta.Z();
  vtx[nPar-1]=new_vtx;
}
