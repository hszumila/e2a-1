/////////////////////////////////////////////////////////////////////
// 
//   File: TMapTools.h   Header File for TMapTools library
//   
//Begin_Html  
/*  TMapTools library provides TOOLS to get information from

 Calibration Maps from mySQL database server. 
 This library is based on TMySQL library in ROOT
 
  Initial Author :  Gagik Gavalian (UNH)
  Creation Date  :  03/08/2000 

*/
//End_Html
/////////////////////////////////////////////////////////////////////

#ifndef __TMAPTOOLS__
#define __TMAPTOOLS__


#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <unistd.h>

// SQL Includes

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
// ROOT includes

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TQuery.h"
#include "TItemValue.h"
#include "TItem.h"

#include "TIntBuffer.h"
#include "TFloatBuffer.h"
#include "TStringBuffer.h"


#ifndef __MAKEDUMMLIB__
#include "MySQLServer.h"
#endif


#define debugQUERYPREVIEW     0x20
#define debugEXCEPTIONSPRINT  0x80

class TMapTools {

  TString  fHOSTNAME;  // Database HOST name
  TString  fDATABASE;  // Database to use
  TString  fUSERNAME;  // Database User name 
  TString  fPASSWORD;  // Database User Password

  Int_t    fDebugMode; // Debugmode controls printouts given from various methods
  TString  fItemValueVersion; // This contains Velue table version Development or Production

//
//
//


 public:


  TMapTools(); // Default Constructor
  TMapTools(const char *dbHOSTNAME, const char *dbDATABASE,const char *dbUSERNAME); // Constructor with initializing Database Parameters
 virtual  ~TMapTools(); // Default Destructor.


  void  SetDBParams(const char *dbHOSTNAME,const char *dbDATABASE, const char *dbUSERNAME); // Set Database Parameters
  void  SetHost(const char *dbHOSTNAME);          // Set Database Server name
  void  SetDatabase(const char *dbDATABASE);      // Set Database name
  void  SetPassword(const char *dbPASSWORD);      // Set Database user Passowrd

//
//
TString GetTableName(const char *dSystem,const char *dSubSystem, const char *dItem);
void    ScanMaps();   // Prints out on Standart output  all Systems (Maps)
void    ScanSystem(const char *tSystem);  // Prints on Standard output SusSystems  of System tSystem
void    PrintItem(const char *tSystem, const char *tSubSystem,const char *tItem, Int_t tMinRun, Int_t tMaxRun);
void    PrintItemTable(TSQLServer *fbConn,TItem tMapItem, Int_t fcRunNum);
//

Int_t   GetSystemIDEx(TSQLServer *conn, const char *system_name);
Int_t   GetSystemID(TSQLServer *conn, const char *system_name,Int_t *pSystem_ID);
Int_t   GetSubSystemIDEx(TSQLServer *conn, int system_id, const char *subsystem_name);
Int_t   GetSubSystemID(TSQLServer *conn, int system_id, const char *subsystem_name, Int_t *pSubSystem_ID);
void    GetItemEx(TSQLServer *conn, int subsystem_id, const char *item_name, TItem *fMapItem);
Int_t   GetItem(TSQLServer *fbConn,const char *tSystem, const char *tSubSystem,const char *tItem,TItem *tMapItem);
void    GetItemValueEx(TSQLServer *fbConn, TItem tMapItem, Int_t fcRunNum, const char *fcAuthor, TItemValue *tItemValue);
Int_t   GetItemValue(TSQLServer *fbConn, TItem tMapItem, Int_t fcRunNum, const char *fcAuthor, TItemValue *tItemValue);
//
//
//

TSQLResult & SelectValues(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no, const char *fcAuthor, const char *fcValueType);

Int_t   *GetMapInt(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,Int_t *fiBuffer_length, const char *fcAuthor = "*"); // 

void   GetMapIntBuff(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,TIntBuffer *fIntBuffer, const char *fcAuthor = "*"); // 

Float_t   *GetMapFloat(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,Int_t *fiBuffer_length, const char *fcAuthor = "*");

void   GetMapFloatBuff(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,TFloatBuffer *fFloatBuffer, const char *fcAuthor = "*");

Double_t   *GetMapDouble(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,Int_t *fiBuffer_length, const char *fcAuthor = "*");

TString   *GetMapChar(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,Int_t *fiBuffer_length, const char *fcAuthor = "*");


//


void    PrintInfo();
void    SetDebugMode(Int_t DbgMode);
void    SetDbVersion(const char *fVersion);

//
//
//

ClassDef(TMapTools,1) // Class of Tools for Calibration Constant Maps

};

#endif
