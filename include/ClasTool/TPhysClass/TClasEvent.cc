// Autogenerated Class (Source File)
// Author : G.Gavalian
// Date   : Tue Jul 21 09:03:23 EDT 2009
//

#include "TClasEvent.h"


ClassImp(TClasEvent)


TClasEvent::TClasEvent (){
  SetName("Generic");
  SetEventBank("EVNT");
  fClasBeamType = "photon";
  SetTarget(2212,0.93827,TVector3(0.,0.,0.),TVector3(0.,0.,0.));
}

TClasEvent::~TClasEvent (){

}
TClasEvent::TClasEvent ( const char *name)
{
  SetName(name);
  SetEventBank("EVNT");
  fClasBeamType = "photon";
  SetTarget(2212,0.93827,TVector3(0.,0.,0.),TVector3(0.,0.,0.));
}
TClasEvent::TClasEvent ( const char *name, const char *partinfo, const char *partBank)
{
  SetName(name);
  SetReaction(partinfo);
  SetEventBank(partBank);
  fClasBeamType = "photon";
  SetExclusiveFlag("");
  SetTarget(2212,0.93827,TVector3(0.,0.,0.),TVector3(0.,0.,0.));
}

TClasEvent::TClasEvent ( const TClasEvent &event)
{
  SetName(event.GetName());
  SetReaction(event.GetReaction().Data());
  SetEventBank(event.GetEventBank().Data());
  SetBeamType(event.GetBeamType().Data());
  SetExclusiveFlag(event.GetExclusiveFlag().Data());
}
//-----------------
//-----------------
//-----------------
void     TClasEvent::SetTarget(Int_t  _pdgId, Double_t _mass, TVector3 _p3, TVector3 _v3)
{
  fClasTarget.SetId(_pdgId);
  fClasTarget.SetPdgMass(_mass);
  fClasTarget.SetVector(TLorentzVector(_p3.X(),_p3.Y(),_p3.Z(),
				       sqrt(_p3.Mag2()+_mass*_mass)));
  fClasTarget.SetVertex(_v3);
  fClasTarget.SetStatus(1);
}
/*---------------------------------------------------*/
void     TClasEvent::SetBeam(Int_t _pdgId, Double_t _mass, TVector3 _p3, TVector3 _v3)
{
  fClasBeam.SetId(_pdgId);
  fClasBeam.SetPdgMass(_mass);
  fClasBeam.SetVector(TLorentzVector(_p3.X(),_p3.Y(),_p3.Z(),
				       sqrt(_p3.Mag2()+_mass*_mass)));
  fClasBeam.SetVertex(_v3);
}
/*---------------------------------------------------*/
void     TClasEvent::Reset()
{
  fExcessCharged = 0;
  fExcessNeutral = 0;
  for(UInt_t i=0;i<fPART_Store.size();i++)
    {
      fPART_Store[i].Reset();
    }
}

Double_t  TClasEvent::Get(const char *_par)
{
  return 0.;
}
/*---------------------------------------------------*/
void     TClasEvent::SetExclusiveFlag(const char *_excl_flag)
{
  fExclusiveFlag = _excl_flag;
  if(fExclusiveFlag.Contains("charged"))
    {
      printf("Exclusive Event flag for charged (+/-) particles is set\n");
    }
  if(fExclusiveFlag.Contains("neutral"))
    {
      printf("Exclusive Event flag for neutral (0/0) particles is set\n");
    }
}
/*---------------------------------------------------*/
//*****************************************************************
// Builds the Event Accosring to the Reaction set for the 
// Builder.
//*****************************************************************
void     TClasEvent::BuildEvent(TVirtualReader &fReader)
{
  Reset();
  fTAGRPhoton.TTAG = 0.;
  fTAGRPhoton.T_id = 0;
  fTAGRPhoton.E_id = 0;

  int nrows = fReader.GetNRows(fRootDSTEventBank.Data());

  //--------------------------------------------------------
  // Loop over all particles in the event
  // Fill the fPART_Store vector with the desired particles
  // The particles that are not in the reaction list
  // will be filled in fExcessCharged and fExcessNeutral
  //--------------------------------------------------------
  //  printf("%s Contains %d rows\n",fRootDSTEventBank.Data(),nrows);

  for(int j=0;j<nrows;j++)
    {
      TVirtualData *ptr_DATA = static_cast<TVirtualData *> 
	(fReader.GetBankRow(fRootDSTEventBank.Data(),j));

      //      printf("%s row %d  id=%d stat=%d\n",fRootDSTEventBank.Data(),j,ptr_DATA->GetId(),ptr_DATA->GetStat());
      if(ptr_DATA->GetStat()<=0&&fRootDSTEventBank.Contains("EVNT")) continue;

      int pid = ptr_DATA->GetId();
      int freeIndex = FindFreeIndex(pid);
      //      printf("row %d  id = %d index = %d\n",j,ptr_DATA->GetId(),freeIndex);
      if(freeIndex<0){
	if(ptr_DATA->GetCharge()==0)
	  {
	    fExcessNeutral++;
	  } else fExcessCharged++;
      } else {
	fPART_Store[freeIndex].SetPartById(pid,ptr_DATA->GetMomVec(),
					    ptr_DATA->GetVertex());
	fPART_Store[freeIndex].SetStatus(1);
      }
    }
  //---------------------------------------------------
  // End of the Loop over the particles
  //---------------------------------------------------

  if(fClasBeamType.Contains("photon"))
    FindPhotonBeamTAGR(fReader);
  
}

/*---------------------------------------------------*/
Int_t    TClasEvent::FindFreeIndex(int id)
{
  for(UInt_t i=0;i<fPartInfoPdgId.size();i++)
    {
      if(id==fPartInfoPdgId[i]&&fPART_Store[i].GetStatus()!=1)
	return i;
    }
  return -1;
}
//*****************************************************************
TLorentzVector  TClasEvent::GetMissing(int indx)
{
  if(indx<0||indx>=fPART_Store.size())
    return TLorentzVector(0.,0.,100.,sqrt(10000.+10000.));
  if(fClasBeam.GetStatus()!=1||fClasTarget.GetStatus()!=1||
     fPART_Store[indx].GetStatus()!=1) 
    return TLorentzVector(0.,0.,100.,sqrt(10000.+10000.));
  
  return (fClasBeam.GetVector() + fClasTarget.GetVector() 
	  - fPART_Store[indx].GetVector());
}
//*****************************************************************
TLorentzVector  TClasEvent::GetInvariant(int indx1, int indx2)
{
  if(indx1<0||indx1>=fPART_Store.size()||
     indx2<0||indx2>=fPART_Store.size())
    return TLorentzVector(0.,0.,100.,sqrt(10000.+10000.));
  return  (fPART_Store[indx1].GetVector()+fPART_Store[indx2].GetVector());
}


void     TClasEvent::FindPhotonBeamTAGR(TVirtualReader &fR)
{
  fClasBeam.Reset();
  THEADERClass  *ptrHEAD = fR.GetHEADER();
  double  StartTime = ptrHEAD->GetSTT();

  int     bestTagrIndex = -1;
  double  bestTimeDiff  = 100.;
  int  nrows = fR.GetNRows("TAGR");
  for(int i=0;i<nrows;i++)
    {
      TTAGRClass  *ptrTAGR = static_cast<TTAGRClass *> (fR.GetBankRow("TAGR",i));
      double  time_diff = ptrTAGR->GetTagRF() - StartTime;
      if(fabs(time_diff)<fabs(bestTimeDiff))
	{
	  bestTimeDiff  = time_diff;
	  bestTagrIndex = i;
	}
    }

  if(bestTagrIndex>=0)
    {
      TTAGRClass *TAGR_ptr = static_cast<TTAGRClass *> 
	(fR.GetBankRow("TAGR",bestTagrIndex));
      	    fClasBeam.SetId(22);
	    fClasBeam.SetTofMass(0.);
	    fClasBeam.SetPdgMass(0.);
	    fClasBeam.SetVector(TLorentzVector(0.,0.,TAGR_ptr->GetEnergy(),
					       TAGR_ptr->GetEnergy()));
	    fClasBeam.SetVertex(TVector3(0.,0.,-200.));
	    fTAGRPhoton.ERG  = TAGR_ptr->ERG;
	    fTAGRPhoton.TTAG = bestTimeDiff;
	    fTAGRPhoton.T_id = TAGR_ptr->T_id;
	    fTAGRPhoton.E_id = TAGR_ptr->E_id;
	    fClasBeam.SetStatus(1);
    } else {
	  cout << "TClasEvent::FindPhotonBeamSEB: ERROR : "
	       << " The index in TGPB[" << bestTagrIndex 
	       << "] does not exist in TAGR[" 
	       << "]" << endl;
  }
}
//*****************************************************************
// Find best phtoton from TGPB Bank
// 
//*****************************************************************
void     TClasEvent::FindPhotonBeamSEB(TVirtualReader &fR)
{
  fClasBeam.Reset();
  Int_t      nrows_TGPB = fR.GetNRows("TGPB");
  Double_t   fClosestTime  = 999.;
  Int_t      fClosestIndex = -1;
  
  //  printf("Calling the FindPhoton Routine rows=%d\n",nrows_TGPB);

  for(int i=0;i<nrows_TGPB;i++)
    {
      TTGPBClass *tgpb_ptr = static_cast<TTGPBClass *> 
	(fR.GetBankRow("TGPB",i));

      if(tgpb_ptr != NULL)
	{
	  if(fabs(tgpb_ptr->GetDt())<fabs(fClosestTime))
	    {
	      fClosestTime  = tgpb_ptr->GetDt();
	      fClosestIndex = tgpb_ptr->GetPointer() - 1;
	    }
	} else {
	cout << "TClasEvent::FindPhotonBeamSEB: ERRROR "
	     << " TGPB pointer for row " << i 
	     <<  " is NULL" << endl;
      }
    }
  //-
  // At this point the fClosestTime contains the Tagger time closest to
  // 0. and fClosestIndex contains pointer to the TAGR bank
  // The following part checks if the hit exists in TAGR, if not 
  // the photon will not be identified, if yes, the TAGR information
  // will be filled with the best photon.
  //-

  Int_t  idx_d = (int) (fClosestIndex/1000);
  Int_t nrows_TAGR = fR.GetNRows("TAGR");
  if(abs(idx_d)>0&&abs(idx_d)<=nrows_TAGR)
    fClosestIndex = abs(idx_d)-1;
  //  printf("Best Time = %f idx = %d %d [%d]\n", fClosestTime,fClosestIndex,idx_d,nrows_TAGR);
  if(fClosestIndex>=0)
    {
      if(fClosestIndex>=nrows_TAGR){
	cout << "TClasEvent::FindPhotonBeamSEB: ERROR : "
	     << " The index in TGPB[" << fClosestIndex 
	     << "] does not exist in TAGR[" << nrows_TAGR
	     << "]" << endl;
      } else {
	TTAGRClass *TAGR_ptr = static_cast<TTAGRClass *> 
	  (fR.GetBankRow("TAGR",fClosestIndex));
	if(TAGR_ptr!=NULL)
	  {
	    fClasBeam.SetId(22);
	    fClasBeam.SetTofMass(0.);
	    fClasBeam.SetPdgMass(0.);
	    fClasBeam.SetVector(TLorentzVector(0.,0.,TAGR_ptr->GetEnergy(),
					       TAGR_ptr->GetEnergy()));
	    fClasBeam.SetVertex(TVector3(0.,0.,-200.));
	    fTAGRPhoton.ERG  = TAGR_ptr->ERG;
	    fTAGRPhoton.TTAG = fClosestTime;
	    fTAGRPhoton.T_id = TAGR_ptr->T_id;
	    fTAGRPhoton.E_id = TAGR_ptr->E_id;
	    fClasBeam.SetStatus(1);
	  } else{
	  cout << "TClasEvent::FindPhotonBeamSEB: ERROR : "
	       << " The index in TGPB[" << fClosestIndex 
	       << "] does not exist in TAGR[" << nrows_TAGR
	       << "]" << endl;
	}
      }
    }
}

//*****************************************************************
// Returns Event Status if 
//  1 - event passed the selection cuts
// -1 - event did not pass the selection cuts
// -2 - event did not have a good photon
// -4 - not selection cuts no good photon
//*****************************************************************
Int_t    TClasEvent::GetEventStatus()
{
  int  status = 1;
  for(UInt_t i=0;i<fPART_Store.size();i++)
    if(fPART_Store[i].GetStatus()!=1) status = -1;

  //  if(fClasBeam.GetStatus()!=1) status -= 3;

  return status;
}
//*****************************************************************
// Set Reaction to be studied
//*****************************************************************
void     TClasEvent::SetReaction(const char *partinfo)
{
  TString str_line = partinfo;
  fPartInfo        = partinfo;

  TObjArray *fList = str_line.Tokenize(":");

  fPART_Store.resize(fList->GetEntries());
  fPartInfoPdgId.resize(fList->GetEntries());
  printf("---\n");
  printf("Event Builder %s Define Reaction\n",GetName());

  for(int i=0;i<fList->GetEntries();i++)
    {
      fPartInfoPdgId[i] = (static_cast<TObjString *> (fList->At(i)))->GetString().Atoi();
      
      printf(" -> Add Particle (id) = %d\n",fPartInfoPdgId[i]);
      fPART_Store[i].Reset();
    }
}

//*****************************************************************
// Print Event information
//*****************************************************************

void     TClasEvent::PrintInfo()
{
  
}

void     TClasEvent::PrintEvent()
{
  printf("CLASEVENT  Run %d  %d (E = %9.5f T_id=%d, E_id=%d, time=%9.5f)\n",
	 fClasRunNumber,fClasEventNumber, fClasBeam.GetVector().P(),
	 fTAGRPhoton.T_id,fTAGRPhoton.E_id, fTAGRPhoton.TTAG);
  printf("Excess neutral/charged  %d/%d\n\n",fExcessNeutral,
	 fExcessCharged);

  for(UInt_t i=0;i<fPART_Store.size();i++)
    {
      printf("%s\n",fPART_Store[i].ToString().Data());
    }
  printf("\n");
}


TClasEvent &TClasEvent::operator=(const TClasEvent &event)
{
  SetName(event.GetName());
  SetReaction(event.GetReaction().Data());
  SetEventBank(event.GetEventBank().Data());
  SetBeamType(event.GetBeamType().Data());
  SetExclusiveFlag(event.GetExclusiveFlag().Data());  
  return *this;
}