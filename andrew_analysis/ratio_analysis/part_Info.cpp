#include "part_Info.h"
part_Info::part_Info(int particleID, double mom_x, double mom_y, double mom_z, double vertex)
{
    parID=particleID;
    mom.SetXYZ(mom_x,mom_y,mom_z);
    vtx=vertex;
}

part_Info::~part_Info()
{
}

//Particle ID functions
bool part_Info::isElectron(){
  if(parID == eCode){
    return true;
  }
  return false;
}


bool part_Info::isNucleon(){
  if(parID == pCode){
    return true;
  }
  else if(parID == nCode){
    return true;
  }
  return false;
}

bool part_Info::isProton(int j){
  int parID = parparID[j];
  if(parID == pCode){
    return true;
  }
  return false;
}

bool part_Info::isNeutron(int j){
  int parID = parparID[j];
  if(parID == nCode){
    return true;
  }
  return false;
}


bool part_Info::isPion(int j){
  int parID = parparID[j];
  if(parID == pipCode){
    return true;
  }
  else if(parID == pimCode){
    return true;
  }
  else if(parID == pi0Code){
    return true;
  }
  return false;
}

bool part_Info::isDelta(int j){
  int parID = parparID[j];
  if(parID == dppCode){
    return true;
  }
  else if(parID == dpCode){
    return true;
  }
  else if(parID == d0Code){
    return true;
  }
  else if(parID == dmCode){
    return true;
  }
  return false;
}

//Functions to get numbers
int part_Info::getParID()
{
  return parID;
}

double part_Info::getPX()
{
  return mom.X();
}

double part_Info::getPY()
{
  return mom.Y();
}

double part_Info::getPZ()
{
  return mom.Z();
}

double part_Info::getVTX()
{
  return vtx[i];
}

TVector3 part_Info::getVector()
{
  return mom;
}
