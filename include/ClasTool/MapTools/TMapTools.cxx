///////////////////////////////////////////
//
//  File :  TMapTools.cxx
//
//  TMapTools is set of tools to get or put
//  data into CLAS Calibration Database (MySQL version)
//
//  Author :   Gagik Gavalian (UNH) 
//  Date   :   04/20/2000
//
///////////////////////////////////////////

#include "TMapTools.h"

ClassImp(TMapTools)


TMapTools::TMapTools(){
  // Default Destructor
  // 

fPASSWORD = "";
SetDbVersion("DEV");

}


TMapTools::TMapTools(const char *dbHOSTNAME, const char *dbDATABASE,const char *dbUSERNAME){
  // 
  // Constructor initializes Database host
  // Database and user name

  fHOSTNAME = dbHOSTNAME;
  fDATABASE = dbDATABASE;
  fUSERNAME = dbUSERNAME;
  fPASSWORD = "";
  SetDbVersion("DEV");

}

TMapTools::~TMapTools(){
// Default destructor
//
//

fDebugMode = 0;

}

void TMapTools::SetHost(const char *dbHOSTNAME){
// Sets mySQL Server host name to dbHOSTNAME
//
fHOSTNAME = dbHOSTNAME;
}

void TMapTools::SetDBParams(const char *dbHOSTNAME,const char *dbDATABASE, const char *dbUSERNAME){
// Sets Database Information i.e.
// HOSTNAME, DATABASE, USERNAME
//
  fHOSTNAME = dbHOSTNAME;
  fUSERNAME = dbUSERNAME;
  fDATABASE = dbDATABASE;

}

void TMapTools::PrintInfo(){
// Prints information about Database settings
// HOSTNAME, DATABASE, USERNAME 
//
//
cout << setiosflags(ios::left);
cout << endl << " Class TMap Tools : Database Information " << endl;
cout << "+-----------------+---------------------------------+"<<endl;
  cout << "| " << setw(16) << " Database  Host " << "| " << setw(32) << fHOSTNAME.Data() << "|" << endl;
cout << "+-----------------+---------------------------------+"<<endl;
  cout << "| " << setw(16) << " DataBase " << "| " << setw(32) << fDATABASE.Data() << "|" << endl;
cout << "+-----------------+---------------------------------+"<<endl;
  cout << "| " << setw(16) << " Database User " << "| " << setw(32) << fUSERNAME.Data() << "|" << endl;
cout << "+-----------------+---------------------------------+"<<endl;
  cout << "| " << setw(16) << " Password " << "| " << setw(32) << fPASSWORD.Data() << "|" << endl;
cout << "+-----------------+---------------------------------+"<<endl<<endl;

}


void  TMapTools::ScanMaps(){
//  Connectos to DATBASE on HOSTNAME and Prints out System Names from
//  table System. They correcspond to MAP names in old version of MapManager
//

 TString serverHost;
 serverHost += "mysql://";
 serverHost += fHOSTNAME;

 TSQLServer *db = TSQLServer::Connect(serverHost.Data(),fUSERNAME.Data(),fPASSWORD.Data());
 if(db->SelectDataBase(fDATABASE.Data()) != 0){
    cout << " ScanMaps: Error selecting DataBase " << fDATABASE.Data() << endl;
    return;
 }

 TSQLRow     *fRow;
 TSQLResult  *fResult;
 TQuery       fQuery;

      cout << setiosflags(ios::left);

  if(db->IsConnected()){


      fQuery << "select * from System";
      fResult = db->Query(fQuery.Data());

	cout << endl << endl << "+--------+-------------------------+" << endl;
	cout << "|  SysID |" << setw(25) << "  System Name  " << "|" << endl;
	cout << "+--------+-------------------------| " << endl;
      while((fRow=fResult->Next())){
	cout << "|  " << setw(4) << fRow->GetField(0) << "  " << setw(1) << "|" << "    " << setw(21) << fRow->GetField(1) << "|" << endl;
      }
	cout << "+--------+-------------------------+" << endl << endl;	
  }

delete db;
}


void TMapTools::ScanSystem(const char *tSystem){
//
// Prints out Subsystems and Items for Given System
//
//

 TString serverHost;
 serverHost += "mysql://";
 serverHost += fHOSTNAME;

 TSQLServer *db = TSQLServer::Connect(serverHost.Data(),fUSERNAME.Data(),fPASSWORD.Data());

 if(db->SelectDataBase(fDATABASE.Data()) != 0){
    cout << " ScanMaps: Error selecting DataBase " << fDATABASE.Data() << endl;
    return;
 }


 TSQLRow     *fRow;
 TSQLRow     *fItemRow; 
 TSQLResult  *fResult;
 TSQLResult  *fItemResult;	
 TQuery       fQuery;
 TQuery       fItemQuery;
 Int_t        tSysID;
// TItemClass   tmpIndex;


  if(db->IsConnected()){
      
      
       GetSystemID(db,tSystem,&tSysID);

      cout << setiosflags(ios::left);

      fQuery << "select subsystemId,subsystemName,description from SubSystem where systemId=" << tSysID;
      fResult = db->Query(fQuery.Data());
      cout << endl << endl << "    ScanMap For System : " << tSystem << endl << endl;
		
      while((fRow=fResult->Next())){
	Int_t  tSubSysID = fQuery.IntOf(fRow->GetField(0));
	cout  << " SubSystem : " << fRow->GetField(1) 
              << " (ID=" << tSubSysID  << ")  Description : " << fRow->GetField(2)  << endl ; 

	fQuery.Reset();
	fQuery << "select * from Item where subsystemId=" << tSubSysID;

	if((fDebugMode&0x20)>0) fQuery.Preview();

	fItemResult = db->Query(fQuery.Data());
	while((fItemRow=fItemResult->Next())){
	  cout << setw(15) << "   " 
               << setw(12) << fItemRow->GetField(1) << " : " 
	       << setw(7) << fItemRow->GetField(0) <<" : " 
	       << setw(7) << fItemRow->GetField(3) << "  :  " 
	       << setw(7) << fItemRow->GetField(4) << "  :  " 
	       << setw(7) << fItemRow->GetField(5) << "  :  " 
	       << endl;
	}
	cout << endl << endl;
      }
  }

delete db;
}


void  TMapTools::PrintItemTable(TSQLServer *fbConn,TItem tMapItem, Int_t fcRunNum){
//
//   
//
//
//
TSQLRow     *fRow;
TSQLResult  *fResult;  
TQuery       fQuery;

cout << endl << "----------------------------------------------------------------------------------------" << endl;
	fQuery << "select * from " << fItemValueVersion.Data() << " where itemId=" << tMapItem.fItemId 
          << " AND " << fcRunNum << ">=minRun AND " <<  fcRunNum << "<=maxRun  ORDER BY time" ;
   fResult = fbConn->Query(fQuery.Data());
   for(int j = 0; j < fResult->GetRowCount(); j++){
	 fRow = fResult->Next();
	 for(int i = 0 ; i < fResult->GetFieldCount(); i++)
	 cout << "  " << fRow->GetField(i);
	 cout << endl;
   }
cout << endl << "----------------------------------------------------------------------------------------" << endl;
}

void  TMapTools::PrintItem(const char *tSystem, const char *tSubSystem,const char *tItem, Int_t tMinRun, Int_t tMaxRun){
// 
//  This function Prints content of ValueDescription table for Given Item
// specified by (System,SubSystem,Item) in Run ranges tMinRun -> tMaxRun
//
// NOTE. This will print all entries independent of author entry in the database
//

 TString serverHost;
 serverHost += "mysql://";
 serverHost += fHOSTNAME;

 TSQLRow     *fRow;
 TSQLResult  *fResult;
 TQuery       fQuery;
 TItem       flMapItem;

 TSQLServer *db = TSQLServer::Connect(serverHost.Data(),fUSERNAME.Data(),fPASSWORD.Data());


if(db->SelectDataBase(fDATABASE.Data()) != 0){
    cout << " PrintTable: Error selecting DataBase " << fDATABASE.Data() << endl;
    return;
 }



 if(GetItem(db,tSystem,tSubSystem,tItem,&flMapItem)==0){

 cout << "--------------------- System found " << endl;
   fQuery << "select * from " << fItemValueVersion.Data() << " where itemId=" << flMapItem.fItemId 
          << " AND minRun>=" << tMinRun << " AND maxRun<=" << tMaxRun;
   	if((fDebugMode&0x20)>0) fQuery.Preview();
   fResult = db->Query(fQuery.Data());
   while((fRow=fResult->Next())){
	cout << setw(4) << "   "
               << setw(12) << fRow->GetField(0) << " : " 
	       << setw(7) << fRow->GetField(2) <<" : " 
	       << setw(7) << fRow->GetField(3) << "  :  " 
	       << setw(7) << fRow->GetField(4) << "  :  " 
	       << setw(7) << fRow->GetField(5) << "  :  " 
	       << setw(7) << fRow->GetField(6) << "  :  " 
	       << endl;		
   }
 }
delete db;
}



TSQLResult & TMapTools::SelectValues(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no, const char *fcAuthor, const char *fcValueType){
//
// Selects Values from Database for given Item , Run Number and author and returns 
// TSQLResult to calling function. Also makes check on type of Item.
//
// NOTE. If calling this function one must call it from "try"-"catch" block.
// If error accures function will trow an exception. 
//

 TQuery  tmpQuery;
 TSQLResult *flResult;
 TItem flMapItem;
 TItemValue flItemValue;

 TString serverHost;
 serverHost += "mysql://";
 serverHost += fHOSTNAME;

 TSQLServer *dbConn = TSQLServer::Connect(serverHost.Data(),fUSERNAME.Data(),fPASSWORD.Data());

 if(dbConn->IsConnected()){

 if(dbConn->SelectDataBase(fDATABASE.Data()) != 0){
    delete dbConn;
	   TQuery localErrMsg;
           localErrMsg << " SelectValues: Error selecting DataBase " << fDATABASE.Data();
	   throw localErrMsg;
 }
 
      if(GetItem(dbConn,fcSystem,fcSubSystem,fcItem,&flMapItem) != 0){
	   delete dbConn;
	   TQuery localErrMsg;
           localErrMsg << "SelectValues: Item Not Found ";
	   throw localErrMsg;
	}

        if(GetItemValue(dbConn,flMapItem,fcRun_no,fcAuthor,&flItemValue) != 0){
	   delete dbConn;
	   TQuery localErrMsg;
           localErrMsg << "SelectValues: Item Not Found ";
	   throw localErrMsg;
	}

             
      if(flMapItem.fType != fcValueType){
	  delete dbConn;
	 TQuery localErrMsg;
         localErrMsg << "SelectValues: Type Is not what was expected   " << fcValueType;
	 throw localErrMsg;
	}


	tmpQuery.Reset(); 
        TQuery  flTblName;
        flTblName << fcSystem << "_" << fcSubSystem << "_" << fcItem << "_VALUE";

	tmpQuery << "select * from " <<  flTblName.Data() << " where itemValueId=" << flItemValue.fItemValueId;
        if((fDebugMode&0x20)>0) tmpQuery.Preview();

        flResult = dbConn->Query(tmpQuery.Data());
        if(flResult->GetRowCount() != 1){
	 delete dbConn;
	 TQuery localErrMsg;
         localErrMsg << "SelectValues: Got unexpected count of Rows " << flResult->GetRowCount();
	 throw localErrMsg;    
        }

	if(flResult->GetFieldCount() != flMapItem.fLength+1){
	 delete dbConn;
	 TQuery localErrMsg;
         localErrMsg << "SelectValues: Wrong Item Length. expected : " << flMapItem.fLength 
                     << "  got : " << flResult->GetFieldCount()-1;
	 throw localErrMsg;
	}       		
} 

return (*flResult);
}

Int_t    *TMapTools::GetMapInt(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no, Int_t *fiBuffer_length, const char *fcAuthor = "*"){
//
// Return Pointer to an Array of Integer Elements extracted from 
// map for given Item. fiBuffer_length return length of array
//
// NOTE. The array is created with new function. user should take care of
//       deleting it after use.
//

Int_t    *flBuffer = NULL; 
TSQLRow  *sqlRow;
TQuery    tmpQuery;

try{

  TSQLResult & sqlResult = SelectValues(fcSystem,fcSubSystem,fcItem,fcRun_no,fcAuthor,"int");
  Int_t lBufferLength = sqlResult.GetFieldCount() -1;
  flBuffer = new Int_t[lBufferLength];
  *fiBuffer_length = lBufferLength;
  sqlRow = sqlResult.Next();
  for(int j=0; j< lBufferLength ; j++) flBuffer[j] = tmpQuery.IntOf(sqlRow->GetField(j+1));

} catch  (TQuery & fErrMsg) {

	cout << endl << fErrMsg.Data() << endl << endl;
	*fiBuffer_length = 0;
	return NULL;

}

return flBuffer;
}

void    TMapTools::GetMapIntBuff(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,TIntBuffer *fIntBuffer, const char *fcAuthor = "*"){
//
//  This routine initializes Integer Buffer passed to it and fills it with  
//  Values for given Item. There are several advantages using this method.
//
//  1) User does not have to warry about deleting array if it is defined
//      as static. For example :
//   {
//       TIntBuffer  myDC_T0;
//       GetMapIntBuff("DC_","t0","value",15678,&myDC_T0);
//   }
//   after this braket buffer will be delated automaticaly.
//
//   2) Class TIntBuffer has bound check (unlike array of int's).
//      if one will ask for non existing element in array it will not 
//      couse segmentation falt.
//
//  NOTE: It is recomended to use this function instead of 
//  GetMapInt to avoid memory leaks. 
// 
// 

TSQLRow  *sqlRow;
TQuery    tmpQuery;

try{

  TSQLResult & sqlResult = SelectValues(fcSystem,fcSubSystem,fcItem,fcRun_no,fcAuthor,"int");
  Int_t lBufferLength = sqlResult.GetFieldCount() -1;
  fIntBuffer->SetSize(lBufferLength);
//  flBuffer = new Int_t[lBufferLength];
//  *fiBuffer_length = lBufferLength;
  sqlRow = sqlResult.Next();
  for(int j=0; j< lBufferLength ; j++) fIntBuffer->Set(j,tmpQuery.IntOf(sqlRow->GetField(j+1)));

} catch  (TQuery & fErrMsg) {

	cout << endl << fErrMsg.Data() << endl << endl;

}

}


Float_t    *TMapTools::GetMapFloat(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no, Int_t *fiBuffer_length, const char *fcAuthor = "*"){
//
// Return Pointer to an Array of Float Elements extracted from 
// map for given Item. fiBuffer_length return length of array
//
// NOTE. The array is created with new function. user should take care of
//       deletin it after use.


Float_t    *flBuffer = NULL; 
TSQLRow  *sqlRow;
TQuery    tmpQuery;

try{

  TSQLResult & sqlResult = SelectValues(fcSystem,fcSubSystem,fcItem,fcRun_no,fcAuthor,"float");
  Int_t lBufferLength = sqlResult.GetFieldCount()-1;
  flBuffer = new Float_t[lBufferLength];
  *fiBuffer_length = lBufferLength;

  sqlRow = sqlResult.Next();

  for(int j=0; j< lBufferLength ; j++){ 
     flBuffer[j] = tmpQuery.FloatOf(sqlRow->GetField(j+1));
  }

} catch  (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	*fiBuffer_length = 0;
	return NULL;
}

return flBuffer;
}


void    TMapTools::GetMapFloatBuff(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no,TFloatBuffer *fFloatBuffer, const char *fcAuthor = "*"){
//
//  This routine initializes Float Buffer passed to it and fills it with  
//  Values for given Item. There are several advantages using this method.
//
//  1) User does not have to warry about deleting array if it is defined
//      as static. For example :
//   {
//       TFloatBuffer  myDC_T0;
//       GetMapFloatBuff("DC_","t0","value",15678,&myDC_T0);
//   }
//   after this braket buffer will be delated automaticaly.
//
//   2) Class TFloatBuffer has bound check (unlike array of int's).
//      if one will ask for non existing element in array it will not 
//      couse segmentation falt.
//
//  NOTE: It is recomended to use this function instead of 
//  GetMapFloat to avoid memory leaks. 
// 
// 



TSQLRow  *sqlRow;
TQuery    tmpQuery;

try{

  TSQLResult & sqlResult = SelectValues(fcSystem,fcSubSystem,fcItem,fcRun_no,fcAuthor,"float");
  Int_t lBufferLength = sqlResult.GetFieldCount()-1;
  fFloatBuffer->SetSize(lBufferLength);
//  flBuffer = new Float_t[lBufferLength];
//  *fiBuffer_length = lBufferLength;

  sqlRow = sqlResult.Next();

  for(int j=0; j< lBufferLength ; j++){ 
     fFloatBuffer->Set(j,tmpQuery.FloatOf(sqlRow->GetField(j+1)));
  }

} catch  (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
}

}

TString    *TMapTools::GetMapChar(const char *fcSystem, const char *fcSubSystem,const char *fcItem, Int_t fcRun_no, Int_t *fiBuffer_length, const char *fcAuthor = "*"){
//
// Return Pointer to an Array of TString Elements extracted from 
// map for given Item. fiBuffer_length return length of array
//
// NOTE. The array is created with new function. user should take care of
//       deletin it after use.
//

TString    *flBuffer = NULL; 
TSQLRow   *sqlRow;
TQuery    tmpQuery;

try{

  TSQLResult & sqlResult = SelectValues(fcSystem,fcSubSystem,fcItem,fcRun_no,fcAuthor,"char");
  Int_t lBufferLength = sqlResult.GetFieldCount()- 1;
  flBuffer = new TString[lBufferLength];
  *fiBuffer_length = lBufferLength;

  sqlRow = sqlResult.Next();

  for(int j=0; j< lBufferLength ; j++){ 
     flBuffer[j] = sqlRow->GetField(j+1);
  }

} catch  (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	*fiBuffer_length = 0;
	return NULL;
}

return flBuffer;
}


Int_t    TMapTools::GetSystemID(TSQLServer *conn, const char *system_name, Int_t *pSystem_ID){
//
// This routine writes into pSystem_ID system ID for given System_name
//
// returns 0 if system was found. -1 otherwise
//

try{
   *pSystem_ID = GetSystemIDEx(conn,system_name);
}catch  (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	return -1;
}
 return 0;
}

Int_t    TMapTools::GetSystemIDEx(TSQLServer *conn, const char *system_name){
Int_t System_id = -1;

   TSQLRow     *fRow;
   TSQLResult  *fResult;
   TQuery fQuery;

	fQuery << "select systemId from System where systemName=\"" << system_name << "\"";
	if((fDebugMode&0x20)>0) fQuery.Preview();
	fResult = conn->Query(fQuery.Data());

	if(fResult->GetRowCount()==1){
	   fRow = fResult->Next();
	   sscanf(fRow->GetField(0),"%d",&System_id);
	} else {

	   TQuery ErrMsg;
           ErrMsg << "==>> ERROR(TItemClass::GetSystemID): System " << system_name << " not Found ";
	   throw ErrMsg;
	}

 return System_id;
}


Int_t    TMapTools::GetSubSystemID(TSQLServer *conn, int system_id, const char *subsystem_name, Int_t *pSubSystem_ID){
//
//  SubSystem in pSubSystem_ID variable
//
//

try{
   *pSubSystem_ID = GetSubSystemIDEx(conn,system_id,subsystem_name);
}catch  (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	return -1;
}
return 0;
}

Int_t    TMapTools::GetSubSystemIDEx(TSQLServer *conn, int system_id, const char *subsystem_name){
//
// Returns SubSystem id from table SubSystem corresponding to given SystemID
//
//

Int_t SubSystem_id = -1;

   TSQLRow     *fRow;
   TSQLResult  *fResult;
   TQuery fQuery;
   
        fQuery << "select subsystemId from SubSystem where systemId=" 
               << system_id << " AND subsystemName=\"" << subsystem_name <<"\"";

	if((fDebugMode&0x20)>0) fQuery.Preview();

	fResult = conn->Query(fQuery.Data());

	if(fResult->GetRowCount()==1){
	   fRow = fResult->Next();
	   sscanf(fRow->GetField(0),"%d",&SubSystem_id);
//	   cout << " System ID = " << fRow->GetField(0) << endl;
	} else {
	   TQuery ErrMsg;
           ErrMsg << "==>> ERROR(TItemClass::GetSubSystemID): SubSystem not found. System ID = " << system_id 
                  << "  SubSystem Name : " << subsystem_name;
	   throw ErrMsg;
	}

return SubSystem_id;
}

void    TMapTools::GetItemEx(TSQLServer *conn, int subsystem_id, const char *item_name, TItem *fMapItem){

   TSQLRow     *fRow;
   TSQLResult  *fResult;
   TQuery       fQuery;

	fQuery << "select * from Item where subsystemID=" 
               << subsystem_id << " AND itemName=\"" << item_name << "\"";
	if((fDebugMode&0x20)>0) fQuery.Preview();

	fResult = conn->Query(fQuery.Data());

	if(fResult->GetRowCount()==1){
	   fRow = fResult->Next();
	   
	   fMapItem->fItemId        = fQuery.IntOf(fRow->GetField(0));
	   fMapItem->fItemName      = (char *) fRow->GetField(1);
	   fMapItem->fSubsystemId   = fQuery.IntOf(fRow->GetField(2));
           fMapItem->fLength        = fQuery.IntOf(fRow->GetField(3));
	   fMapItem->fType          = (char *) fRow->GetField(4);
	   fMapItem->fDescription   = (char *) fRow->GetField(5);

	} else {
	   TQuery ErrMsg;
           ErrMsg << "==>> ERROR(TItemClass::GetItemID): Item not found. SubSystem ID = " << subsystem_id 
                  << "  Item Name : " << item_name;
	   throw ErrMsg;
	}
}


Int_t  TMapTools::GetItem(TSQLServer *fbConn, const char *fcSystem,const char *fcSubSystem,const char *fcItem, TItem *sfItem){

 try{

	Int_t  tmpSystemID    = GetSystemIDEx(fbConn,fcSystem);
        Int_t  tmpSubSystemID = GetSubSystemIDEx(fbConn,tmpSystemID,fcSubSystem);
	GetItemEx(fbConn,tmpSubSystemID,fcItem,sfItem);

 } catch (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	return -1;
 } 

return 0;
}

Int_t   TMapTools::GetItemValue(TSQLServer *fbConn, TItem tMapItem, Int_t fcRunNum, const char *fcAuthor, TItemValue *tItemValue){
try{
     GetItemValueEx(fbConn,tMapItem,fcRunNum,fcAuthor,tItemValue);
} catch (TQuery & fErrMsg) {
	cout << endl << fErrMsg.Data() << endl << endl;
	return -1;
}
return 0;

}
void   TMapTools::GetItemValueEx(TSQLServer *fbConn, TItem tMapItem, Int_t fcRunNum, const char *fcAuthor, TItemValue *tItemValue){

   TSQLRow     *fRow;
   TSQLResult  *fResult;
   TQuery       fQuery;

   fQuery.Reset();
if(strstr(fcAuthor,"*") != NULL){
	cout << "Getting Thigs Ordered by time " << endl;
	fQuery << "select * from " << fItemValueVersion.Data() << " where itemId=" << tMapItem.fItemId 
          << " AND " << fcRunNum << ">=minRun AND " <<  fcRunNum << "<=maxRun  ORDER BY time" ;
} else {
   
   fQuery << "select * from " << fItemValueVersion.Data() << " where itemId=" << tMapItem.fItemId
          << " AND author=\""<< fcAuthor <<"\" AND " << fcRunNum 
          << ">=minRun AND " <<  fcRunNum << "<=maxRun ORDER BY time";
}

  if((fDebugMode&0x20)>0) fQuery.Preview();

   fResult = fbConn->Query(fQuery.Data());

   if(fResult->GetRowCount()<=0){
	 cout << "GetMapValues: No record where found. Possible candidates are: " << endl;
	 PrintItemTable(fbConn,tMapItem,fcRunNum);
	 TQuery ErrMsg;
           ErrMsg << "==>> ERROR(TMapTools::GetItemValueID): Item Value not found. Item ID = " << tMapItem.fItemId 
                  << " Run : " << fcRunNum << "  Author : " << fcAuthor;
	   throw ErrMsg;   
   } else {
	   if(fResult->GetRowCount() > 1){
		cout << endl << endl << "\t\t ==>> W A R N I N G <<== " << endl << endl
                     << "GetItemValues: " << fResult->GetRowCount() << " records were found. Using most recent one. "
                     << endl << endl;
	}
	   for(int j = 0; j < fResult->GetRowCount(); j++) fRow = fResult->Next();
//	   fRow = fResult->Next();
	   tItemValue->fItemValueId = fQuery.IntOf(fRow->GetField(0));
	   tItemValue->fItemId  = fQuery.IntOf(fRow->GetField(1));
	   tItemValue->fMinRun  = fQuery.IntOf(fRow->GetField(2));
	   tItemValue->fMaxRun  = fQuery.IntOf(fRow->GetField(3));
	   tItemValue->fAuthor  = fRow->GetField(4);
	   tItemValue->fTime    = fRow->GetField(4);
	   tItemValue->fComment = fRow->GetField(4);
   }
	  
}


TString TMapTools::GetTableName(const char *dSystem,const char *dSubSystem, const char *dItem){

TString  tblName;

tblName = dSystem;
tblName += "_" ;
tblName = dSubSystem;
tblName += "_" ;
tblName = dItem;
tblName += "_VALUE" ;

return tblName;
}


void TMapTools::SetDebugMode(Int_t DbgMode){
fDebugMode = DbgMode;
}


void    TMapTools::SetDbVersion(const char *fVersion){
//
//
//
//
  if(strcmp(fVersion,"DEV") == 0){
	fItemValueVersion = "ItemValueDescription";
	cout << "SetGbVersion: Set Database version to Development. Table \"" 
             << fItemValueVersion.Data() <<"\" will be used" << endl;
  }
  if(strcmp(fVersion,"PROD") == 0){
	fItemValueVersion = "ProductionItemValueDescription";
	cout << "SetGbVersion: Set Database version to Production. Table \"" 
             << fItemValueVersion.Data() <<"\" will be used" << endl;
  }
}
