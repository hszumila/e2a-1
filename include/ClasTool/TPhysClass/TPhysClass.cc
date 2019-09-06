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


#include "TPhysClass.h"

ClassImp(TPhysClass)

TPhysClass::TPhysClass(){
  cout << "TPHYSCLASS: Initializing the class" << endl;
}

TPhysClass::~TPhysClass()
{
}

void  TPhysClass::AddEventBuilder(TClasEvent *evtBuilder)
{

  TObject *check_obj = fEventBuilderStore.FindObject(evtBuilder->GetName());
  
  if(check_obj!=NULL)
    {
      cout << "AddEventBuilder: Failed to Add TClasEvent "
	   << "with name \"" << evtBuilder->GetName() 
	   << "\". The name already exists" << endl;
      return;
    }
  fEventBuilderStore.Add(evtBuilder);
}

int  TPhysClass::NextEvent(const char *bank)
{

  Int_t  status = Next();
  int nrow = GetNRows(bank);

  for(int kk=0;kk<fEventBuilderStore.GetEntries();kk++)
    {
      TClasEvent *evtBuilder = static_cast<TClasEvent *>
	(fEventBuilderStore.At(kk));
      evtBuilder->BuildEvent(*this);
    }

  return status;
}

TClasEvent   *TPhysClass::GetEventBuilder(const char *evtbName)
{
  return static_cast<TClasEvent *> (fEventBuilderStore.FindObject(evtbName));
}

void   TPhysClass::PrintEventBank(const char *bankname )
{
  int nrows = GetNRows(bankname);
  printf("---  Printout Bank %s - Run %d  %d  --\n",bankname,
	 GetRunNumber(), GetEventNumber());
  for(int i=0;i<nrows;i++)
    {
      TVirtualData *data = static_cast<TVirtualData *>( GetBankRow(bankname,i));
      if(data->GetStat()>=0)
	printf("%2d %5d  %9.5f  %9.5f  %9.5f  %5d\n",i,
	       data->GetId(),data->GetMomVec().X(),
	       data->GetMomVec().Y(),data->GetMomVec().Z(),data->GetStat());
    }
}
