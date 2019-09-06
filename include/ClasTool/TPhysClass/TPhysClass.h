//=================================================
//
//  TPhysClass is a template class for introducing
//  our own analyses in TClasTool framework
//  just follow examples in this directory
//
//   Author: G.Gavalian (UNH) 04/12/2002
//   
//
//==================================================


#ifndef __TPHYS_CLASS__
#define __TPHYS_CLASS__

#include <iostream>
#include <TObject.h>
#include <TDSTReader.h>
#include <TPartCounter.h>
#include <TObjArray.h>
#include <TClasEvent.h>

using namespace std;

class TPhysClass : public TDSTReader {

  TObjArray     fEventBuilderStore;

 public:

  TPhysClass();
  virtual ~TPhysClass();

  void   AddEventBuilder(TClasEvent *evtBuilder);
  int    NextEvent(const char *bank = "EVNT");
  TClasEvent   *GetEventBuilder(const char *evtbName);
  void   PrintEventBank(const char *bankname = "EVNT");


  ClassDef(TPhysClass,2)
};

#endif
