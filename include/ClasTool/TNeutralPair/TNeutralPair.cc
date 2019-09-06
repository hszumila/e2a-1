////////////////////////////////////////////////////////////////////////
// File TNeutralPair.cc
//
// Initial Author: Hovanes Egiyan  UNH  03/24/06
//
// Class to find best Neutral pair
// This class was written with EG3 data set in mind
// In particular there are a lot of things which 
// are particular for photon runs and will not work 
// for electron runs.
// 
////////////////////////////////////////////////////////////////////////


#include "TNeutralPair.h"

Double_t TNeutralPair::npSigmaMass = 0.0022;      // A little tight
Double_t TNeutralPair::npSigmaDOCA = 1.0 ;        // Approximately right
Double_t TNeutralPair::npSigmaTime = 0.300 ;      // Checked with data

Double_t TNeutralPair::npMassCut = 2.5;           // Number of sigmas for Mass cut
Double_t TNeutralPair::npDOCACut = 2.5;           // Number of sigmas for DOCA cut
Double_t TNeutralPair::npTimeCut = 2.5;           // Number of sigmas for the dT cut

Double_t TNeutralPair::npTagrRes = 0.2;           // Tagger time resolution in ns

// Get values for some particles often used here

Double_t TNeutralPair::npClight = 29.979;        // Speed of light in vacuum

Double_t TNeutralPair::npMinMomPos = 0.200 ;     // Minimum momentum for a track to have
Double_t TNeutralPair::npMinMomNeg = 0.100 ;     // Minimum momentum for a track to have

Int_t   TNeutralPair::npElossRunNo = -1;        // Run number of last pair with eloss
                                                // needed for eloss initialization 

Double_t TNeutralPair::npTargetCenter = -50.0;  // Center of the target in the CLAS during the run.
Double_t TNeutralPair::npTagrErgShift = 0.012;  // Shift in the tagger energy from P. Mattioni

//Double_t TNeutralPair::npTagTimeOff = -0.212;   // Time offset between pair time and RF time
Double_t TNeutralPair::npTagTimeOff = -0.0;   // Time offset between pair time and RF time

Double_t TNeutralPair::npSimLbdTimeOff = 0;   // Time offset between negative and positive tracks in GSIM pair (ns)

TDatabasePDG TNeutralPair::npDbPDG;

TMomCorEG3   TNeutralPair::npMomCor;

TNeutralPair::FitGlobType TNeutralPair::npFitGlob;


// class for tagger info

TNeutralPair::TTagInfo   TNeutralPair::TTagInfo::operator=( TNeutralPair::TTagInfo tag ) {
  Status = tag.Status;
  Sort   = tag.Sort;
  Eid    = tag.Eid;
  Tid    = tag.Tid;
  Energy = tag.Energy;
  Time   = tag.Time;
  return (*this);
}
// End of tagger info class


TNeutralPair::TNeutralPair( TDSTReader* reader,  Bool_t eLossFlag, Bool_t fitFlag, 
			    Int_t ntrId, Int_t posId, Int_t negId, Double_t  mPeak ) :
  npOK(1), npKinFitFlag(0), npElossFlag(0), 
  npPosId( posId ), npNegId( negId ), npNtrId( ntrId ), 
  npXi2(-1), npXi2p(-1), npXi2a(-1), 
  npReader( reader ),  npNeutralPeak( mPeak ), npTag()
{
  // This is the pair CTOR

  npMassPos = npDbPDG.GetParticle( npPosId )->Mass() ;
  npMassNeg = npDbPDG.GetParticle( npNegId )->Mass() ;
  npMassNtr = npDbPDG.GetParticle( npNtrId )->Mass() ;

  if( npNeutralPeak < -1000 ) npNeutralPeak = npMassNtr;

  if ( reader != 0 ) {
    npBeta = new TVectorD( npReader->GetNPart() );
    // If not particles in the event pair is bad 
    if ( npBeta->GetNoElements() < 1 ) 
      npOK = 0;
    
    if ( npOK && FindBestNeutral() > 0 ) { 
      FindBetas();
      if ( eLossFlag ) {
	InitEloss( npReader->GetRunNumber() );
       	if( npReader->GetRunNumber() > 20 ) CorrectMomentum();
	npElossFlag = 1;
      }
      if ( fitFlag ) {
	if ( npFitGlob.MinuitFlag == 0 ) {
	  InitFit(); 
	  npFitGlob.MinuitFlag = 1 ;
	}
	KinFit();
	npKinFitFlag = 1;
      }
    } else {
      // No lambdas found in the event, pair becomes bad
      npOK = 0 ;
    }
  } else {
    // Something wrong with the DST reader, pair becomes bad
    npOK = 0;
    npReader = 0;
    npBeta = 0;
  }
}



TNeutralPair::~TNeutralPair() {
  // DTOR
  if ( npBeta != 0 ) 
    delete npBeta ;
}




TNeutralPair::TNeutralPair( TNeutralPair& Pair ) {
  // Copy CTOR
  npKinFitFlag = Pair.npKinFitFlag ;
  npOK = Pair.npOK;
  npBestPosIdx = Pair.npBestPosIdx ;
  npBestNegIdx = Pair.npBestNegIdx ;
  
  npBestMass         = Pair.npBestMass ;
  npBestDOCA         = Pair.npBestDOCA ;
  npBestVtxTimeDiff  = Pair.npBestVtxTimeDiff ;

  npXi2  = Pair.npXi2;
  npXi2p = Pair.npXi2p;
  npXi2a = Pair.npXi2a;

  npReader = Pair.npReader;  
  npTag    = Pair.npTag;

  if( npReader != 0 ) 
    npBeta = new TVectorD( *(Pair.npBeta) );
  else 
    npBeta = 0;

  npFitMom = Pair.npFitMom;
  npMom = Pair.npMom ;
  npVtx = Pair.npVtx ;  
}



TNeutralPair& TNeutralPair::operator=( TNeutralPair& Pair ) { 
  // Copy operator
  npKinFitFlag = Pair.npKinFitFlag ;
  npOK = Pair.npOK;
  npBestPosIdx = Pair.npBestPosIdx ;
  npBestNegIdx = Pair.npBestNegIdx ;
  
  npBestMass         = Pair.npBestMass ;
  npBestDOCA         = Pair.npBestDOCA ;
  npBestVtxTimeDiff  = Pair.npBestVtxTimeDiff ;

  npXi2  = Pair.npXi2;
  npXi2p = Pair.npXi2p;
  npXi2a = Pair.npXi2a;

  npReader = Pair.npReader;  
  (*npBeta) = (*(Pair.npBeta));

  npTag = Pair.npTag;

  npFitMom = Pair.npFitMom;
  npMom = Pair.npMom ;
  npVtx = Pair.npVtx ;

  return static_cast<TNeutralPair&>(*this) ;
}



Int_t TNeutralPair::FindBestNeutral() { 
  //
  // Loop through all q+q- combinations and 
  // find the best candidate for Neutral
  // Momentum correction are not supposed to be 
  // applied at this point. 
  //
  Int_t iPair = 0;
  if ( npReader != 0 && npOK ) {
    THEADERClass* header = npReader->GetHEADER();

    Int_t N_trk = npReader->GetNPart() ;
    
    Double_t DiffMin = 1.0e+30;

    // Loop for positive candidate first
    for ( Int_t iTrkPos= 0; iTrkPos < N_trk; iTrkPos++ ) {
      // Get the EVNT bank 
      TEVNTClass* evntPos = static_cast<TEVNTClass*>( npReader->GetBankRow( "EVNT", iTrkPos ) );
      if ( evntPos->GetStat() > 0 && evntPos->GetDCStat() > 0 && evntPos->GetSCStat() > 0 && 
	   evntPos->GetMomentum() > npMinMomPos && evntPos->GetCharge() > 0 ) {
	// Get DCPB bank and check if track passed time-based tracking
	TDCPBClass* dcpbPos = static_cast<TDCPBClass*>( npReader->GetBankRow( "DCPB", evntPos->GetDCidx() ) );
	if ( dcpbPos->GetStatus() > 0 ) {
	  // Found positive candidate with TBT
	  TSCPBClass* scpbPos = static_cast<TSCPBClass*>( npReader->GetBankRow( "SCPB", evntPos->GetSCidx() ) );

	  // Fill positive momentum 4-vector
	  TVector3 vP_pos = evntPos->GetMomVec();
	  Double_t dE_pos = sqrt( vP_pos.Mag2() + npMassPos*npMassPos  ) ;
	  TLorentzVector lvPosMom( vP_pos, dE_pos );

	  // Fil positive vertex 4-vector
	  TVector3 vV_pos = evntPos->GetVertex();
	  Double_t dVtxTime = scpbPos->Time - ( scpbPos->Path / ( npClight * lvPosMom.Beta() ) ) ;
	  TLorentzVector lvPosVtx( vV_pos, dVtxTime );

	  // Now loop for negative candidate
	  for ( Int_t iTrkNeg= 0; iTrkNeg < N_trk; iTrkNeg++ ) {
	    // Get the EVNT bank 
	    TEVNTClass* evntNeg = static_cast<TEVNTClass*> ( npReader->GetBankRow( "EVNT", iTrkNeg ) );
	    if ( evntNeg->GetStat() > 0 && evntNeg->GetDCStat() > 0 && evntNeg->GetSCStat() > 0 && 
		 evntNeg->GetMomentum() > npMinMomNeg && evntNeg->GetCharge() < 0 ) {
	      // Get DCPB bank and check if track passed time-based tracking
	      TDCPBClass* dcpbNeg = static_cast<TDCPBClass*>( npReader->GetBankRow( "DCPB", evntNeg->GetDCidx() ) );
	      if ( dcpbNeg->GetStatus() > 0 ) {
		// Found negative candidate with TBT
		iPair++;

		TSCPBClass* scpbNeg = static_cast<TSCPBClass*>( npReader->GetBankRow( "SCPB", evntNeg->GetSCidx() ) );
		
		// Fill negative momentum 4-vector
		TVector3 vP_neg = evntNeg->GetMomVec();
		Double_t dE_neg = sqrt( vP_neg.Mag2() + npMassNeg*npMassNeg  ) ;
		TLorentzVector lvNegMom( vP_neg, dE_neg );
		
		// Fill pion vertex 4-vector
		TVector3 vV_neg = evntNeg->GetVertex();
		Double_t dVtxTime = scpbNeg->Time - ( scpbNeg->Path / ( npClight * lvNegMom.Beta() ) );
		TLorentzVector lvNegVtx( vV_neg, dVtxTime );
		
		// calcualate the lambda candidate 4-mom vector and the mass
		TLorentzVector lvNtrMom = ( lvPosMom + lvNegMom ) ;
		Double_t M_ntr = lvNtrMom.M();
		
		// Calculate neutral candidates vertex 4-vector using 
		// Mauriks DOCA function and pion time 
		Double_t VtxDist = 0;
		TVector3 vV_ntr = FindNeutralVtx( vV_pos, vP_pos, vV_neg, vP_neg, &VtxDist );
		Double_t TimeDiff = lvNegVtx.T() - lvPosVtx.T() ;
		if( npReader->GetRunNumber() < 20 ) TimeDiff -= npSimLbdTimeOff; 


		//		Double_t TimeAvg  = ( lvPimVtx.T() + lvPtrVtx.T() ) / 2.0 ; 	      
		TLorentzVector lvNtrVtx( vV_ntr, lvNegVtx.T() );
		
		// This is an adhoc criteria function for best 
		// Neutral selection 
		Double_t Diff = pow( TMath::Abs( M_ntr - npNeutralPeak ) / npSigmaMass , 2 ) 
		  + pow( VtxDist / (1.5*npSigmaDOCA) , 2 ) + 
		  + pow( TimeDiff / (1.5*npSigmaTime), 2 );
		
		if( Diff < DiffMin ) {
		  // If this pair is better remember the numbers and vectors 
		  // for this pair. 
		  npBestPosIdx = iTrkPos;
		  npBestNegIdx = iTrkNeg;
		  npBestMass = M_ntr;
		  npBestDOCA = VtxDist;
		  npBestVtxTimeDiff = TimeDiff;

		  npNeg4Mom = lvNegMom ;
		  npPos4Mom = lvPosMom ;
		  npMom = lvNtrMom;
		  npVtx = lvNtrVtx;
		  DiffMin = Diff ;
		}
	      }
	    }
	  }
	}
      } 
    }      
    
  } else {
    cerr << "Reader is bad " << endl;
    npOK = 0;
    return 0;
  }

  // If the the minimum cuts are not satisfied for the best lambda candidate then 
  // we assume that there was no lambda found in this event.
  if (  ( TMath::Abs( npBestMass - npNeutralPeak ) > ( npMassCut * npSigmaMass ) ) || 
	TMath::Abs( npBestDOCA ) > ( npDOCACut * npSigmaDOCA ) || 
	TMath::Abs( npBestVtxTimeDiff ) > ( npTimeCut * npSigmaTime ) ) 
    iPair = 0;

  return iPair;
}


Double_t TNeutralPair::FindFlightTime() {
  //
  // Determine flight time of Neutral before decaying
  //
  //
  Double_t Tavg = 0.0 ;
  if ( npReader != 0 && npOK ) {
    Int_t N_trk = npReader->GetNPart() ;
    vector<Double_t> aFltLs; 
    for ( Int_t iTrk = 0; iTrk < N_trk; iTrk++ ) {
      if( iTrk != npBestPosIdx && iTrk != npBestNegIdx ) {
	// Get the EVNT bank 
	TEVNTClass* evnt = static_cast<TEVNTClass*>( npReader->GetBankRow( "EVNT", iTrk ) );
	if ( evnt->GetStat() > 0 && evnt->GetDCStat() > 0 &&  evnt->GetSCStat() > 0  ) {
	  // Get DCPB bank and check if track passed time-based tracking
	  TDCPBClass* dcpb = static_cast<TDCPBClass*>( npReader->GetBankRow( "DCPB", evnt->GetDCidx() ) );
	  if ( dcpb->GetStatus() > 0 && 
	       ( this->GetVtx().Vect().Z() - evnt->GetVertex().Z() ) > -2.0 ) {
	    Double_t dist = ( this->GetVtx().Vect() - evnt->GetVertex() ).Mag() ;
	    Double_t time = dist / ( npMom.Beta() * npClight );
	    aFltLs.push_back( time );
	  } 
	}
      }
    }

    Int_t nGoodTrk = aFltLs.size();
    if ( nGoodTrk > 0 ) {
      for( UInt_t iTrk = 0; iTrk < aFltLs.size(); iTrk++ ) {
	Tavg +=  aFltLs[iTrk] / nGoodTrk ;
      }
    } 
  }
  return Tavg;
}


Double_t TNeutralPair::FindStartTime() {
  if ( npReader != 0 && npOK ) {
    THEADERClass* header = npReader->GetHEADER();
    if( npReader->GetRunNumber() > 20 ) {
      npStartTime = header->GetSTT() -
	TMath::Nint( ( header->GetSTT() - this->GetVtx().T() ) / 2.004 ) * 2.004; 
    } else {
      npStartTime = this->GetVtx().T();
      //      npStartTime = header->GetSTT();
    }
    return npStartTime;
  } else {
    return -9999;
  }
}


void TNeutralPair::FindBetas() {
  //
  // Determine betas for all tracks using the SEB start time 
  // corrected for needed number of 2.004 ns determined from lambda
  //
  //
  if ( npReader != 0 && npOK ) {
    FindStartTime();
    
    //    THEADERClass* header = npReader->GetHEADER();   
 
    // This is the adjusted start time 
    //    npStartTime =  header->GetSTT() - 
    //      TMath::Nint( ( header->GetSTT() - this->GetVtx().T() ) / 2.004 ) * 2.004;

    Int_t N_trk = npReader->GetNPart() ;
    for ( Int_t iTrk= 0; iTrk < N_trk; iTrk++ ) {
      // Get the EVNT bank 
      TEVNTClass* evnt = static_cast<TEVNTClass*>( npReader->GetBankRow( "EVNT", iTrk ) );
      if ( evnt->GetStat() > 0 && evnt->GetDCStat() > 0 &&  evnt->GetSCStat() > 0  ) {
	// Get DCPB bank and check if track passed time-based tracking
	TDCPBClass* dcpb = static_cast<TDCPBClass*>( npReader->GetBankRow( "DCPB", evnt->GetDCidx() ) );
	if ( dcpb->GetStatus() > 0 ) {
	  TSCPBClass* scpb = static_cast<TSCPBClass*>( npReader->GetBankRow( "SCPB", evnt->GetSCidx() ) );
	  
	  // Fill beta arrays with the new numbers
	  (*npBeta)[iTrk] =   scpb->Path / ( ( scpb->Time  - npStartTime ) * npClight );
	} else {
	  // No TBT, use SEB numbers
	  (*npBeta)[iTrk] = evnt->GetBeta();
	}
      } else {
	// Not good event status etc, use SEB beta
	(*npBeta)[iTrk] = evnt->GetBeta();
      }
    }          
  }
}



Double_t TNeutralPair::GetMatchTime( TVector3 pPart, TVector3 vPart ) {
  // This function finds the start time for each particle by using Lambda vertex and 
  // particle vertex. 
  if ( npOK && npReader != 0 ) {
    THEADERClass* header = npReader->GetHEADER();
    Double_t startTime = header->GetSTT();

    TVector3 ntrMom = npMom.Vect();
    TVector3 ntrVtx = npVtx.Vect();
    if( npKinFitFlag ) ntrMom = npFitMom.Vect() ;

    // Fine the common vertex Lambda 
    Double_t dist = -1;
    TVector3 vtxCom = FindNeutralVtx( ntrVtx, ntrMom, vPart, pPart, &dist );
    
    if( dist < ( 3.0 * npSigmaDOCA ) ) {
      Double_t pathNtr = ( ntrVtx - vtxCom ).Mag();
      Double_t betaNtr = ntrMom.Mag() / sqrt( pow( ntrMom.Mag(), 2 ) + npMassNtr*npMassNtr ) ; 
      Double_t timeNtr = this->GetVtx().T() - pathNtr / ( betaNtr * npClight );

      Double_t T_best = 10;
      Int_t idxBest = -1;
      // Loop to find the best photon comparing the photon RF time at the target center 
      // corrected for the neutrol z-position to the neutral vertex time
      Int_t N_TAGR = npReader->GetNRows( "TAGR" );
      for( Int_t idxTAGR = 0; idxTAGR < N_TAGR; idxTAGR++ ) {
	TTAGRClass* tagr = static_cast<TTAGRClass*>( npReader->GetBankRow( "TAGR", idxTAGR ) );
	if( tagr != 0 && ( tagr->GetStat() == 15 || tagr->GetStat() == 7 ) ) {
	  Double_t T_vert = tagr->GetTagRF() + ( vtxCom.Z() - npTargetCenter ) / npClight; 
	  Double_t Diff = TMath::Abs( timeNtr - T_vert ) ;
	  if( Diff < T_best ) {
	    T_best = Diff;
	    idxBest = idxTAGR;
	    startTime = T_vert;
	  }
	}
      }
      // Select the start time if it finds a good photon, otherwise 
      //  use the SEB start time
      if( idxBest < 0 || T_best > 2 ) startTime = header->GetSTT();
    }
    
    return startTime;
    
  } else {
    // Bad pair object
    return 0;
  }
}



TVector3 TNeutralPair::FindNeutralVtx( TVector3 v1, TVector3 p1, 
			    TVector3 v2, TVector3 p2, Double_t *dist ) {
  //
  // Find the mid point on the line segment that connects the
  // two lines (tracks) at the point of closest approach.
  //
  //  Modified funtion from Maurik
  //
  Double_t R1 = 0;
  TVector3 diff = (v1) - (v2);

  // Phat1, Phat2 are the normalized P vectors, which
  // is just the cosine terms.
  // R=Phat1 dot Phat2
  TVector3 p1unit = p1.Unit();
  TVector3 p2unit = p2.Unit();

  Double_t R = p1unit.Dot(p2unit);

  if( R >= 1. ) {
    cout << "Parallel " << R << endl;
    return(-1);
  } else {
    R1 = 1 / ( 1- R*R );
  }

  TVector3 RP21 = p2unit*R - p1unit;
  TVector3 RP12 = p2unit - p1unit*R;

  Double_t DOT21 = diff.Dot(RP21);
  Double_t DOT12 = diff.Dot(RP12);
//
// The M1 (2) vectors point to the closest point on
// track 1 (2) to track 2 (1)
// R1*DOT21 is equal to "alpha"
// R1*DOT12 is equal to "beta"
//

  TVector3 M1 = v1 + p1unit*R1*DOT21;
  TVector3 M2 = v2 + p2unit*R1*DOT12;

  TVector3 diff1 = M2 - M1;
  (*dist) = diff1.Mag() ;

  TVector3 MidPoint = ( M1 + M2 ) * 0.5;

  return  MidPoint;
}



Int_t TNeutralPair::GetTimeBin( Int_t minEcnt, Int_t maxEcnt ) {
  if ( npReader != 0 && npOK ) { 
    if( npTag[minEcnt][maxEcnt].Sort == 0 ) {
      TagTimeSort( minEcnt, maxEcnt );
    } 
    if( npTag[minEcnt][maxEcnt].Sort == 1 && npTag[minEcnt][maxEcnt].Status == 1 ) {
      return npTag[minEcnt][maxEcnt].Tid;
    } 
  }
  return -1;
}


Int_t TNeutralPair::GetErgBin( Int_t minEcnt, Int_t maxEcnt ) {
  if ( npReader != 0 && npOK ) { 
    if( npTag[minEcnt][maxEcnt].Sort == 0 ) {
      TagTimeSort( minEcnt, maxEcnt );
    } 
    if( npTag[minEcnt][maxEcnt].Sort == 1 && npTag[minEcnt][maxEcnt].Status == 1 ) {
      return npTag[minEcnt][maxEcnt].Eid;
    } 
  }
  return -1;
}


Double_t TNeutralPair::GetTagTime( Int_t minEcnt, Int_t maxEcnt ) {
  if ( npReader != 0 && npOK ) { 
    if( npTag[minEcnt][maxEcnt].Sort == 0 ) {
      TagTimeSort( minEcnt, maxEcnt );
    } 
    if( npTag[minEcnt][maxEcnt].Sort == 1 && npTag[minEcnt][maxEcnt].Status == 1 ) {
      return npTag[minEcnt][maxEcnt].Time;
    } 
  }
  return -99999;
}


Double_t TNeutralPair::GetTagEnergy( Int_t minEcnt, Int_t maxEcnt ) {
  if ( npReader != 0 && npOK ) { 
    if( npTag[minEcnt][maxEcnt].Sort == 0 ) {
      TagTimeSort( minEcnt, maxEcnt );
    } 
    if( npTag[minEcnt][maxEcnt].Sort == 1 && npTag[minEcnt][maxEcnt].Status == 1 ) {
      return npTag[minEcnt][maxEcnt].Energy;
    } 
  }
  return 0;
}



Int_t TNeutralPair::TagTimeSort( Int_t minEcnt, Int_t maxEcnt ) {
   if ( npReader != 0 && npOK ) {
    Double_t minDiff = 100;
    Int_t idxBest = -1; 
    //
    // Loop to find the best photon comparing the photon RF time at the target center 
    // corrected for the neutrol z-position to the neutral vertex time
    //
    Int_t N_TAGR = npReader->GetNRows( "TAGR" );
    for( Int_t idxTAGR = 0; idxTAGR < N_TAGR; idxTAGR++ ) {
      TTAGRClass* tagr = static_cast<TTAGRClass*>( npReader->GetBankRow( "TAGR", idxTAGR ) );
      if( tagr != 0 && ( tagr->GetStat() == 15 || tagr->GetStat() == 7 ) && 
	  ( minEcnt <= tagr->GetEid() && tagr->GetEid() <= maxEcnt ) && 
	  this->TidIsOK( tagr->GetTid() ) ) {
 	Double_t T_vert = tagr->GetTagRF() + ( this->GetVtx().Z() - npTargetCenter ) / npClight; 
	//	cout << "Start time is " << npStartTime << " lambda vertex time is " << T_vert << endl;
	Double_t Diff = TMath::Abs( npStartTime - T_vert ) ;
	if( Diff < minDiff ) {
	  minDiff = Diff;
	  idxBest = idxTAGR;
	} 
      } 
    }

    if( idxBest >= 0 && minDiff < 4*npTagrRes ) {
      TTAGRClass* tagr = static_cast<TTAGRClass*>( npReader->GetBankRow( "TAGR", idxBest ) );
      //      npTag[minEcnt][maxEcnt].Energy = tagr->GetEnergy() + npTagrErgShift ;
      npTag[minEcnt][maxEcnt].Energy = TNeutralPair::npMomCor.CorrectEgamma( tagr->GetEnergy() )  ;
      npTag[minEcnt][maxEcnt].Time = tagr->GetTagTime();
      npTag[minEcnt][maxEcnt].Eid = tagr->GetEid();
      npTag[minEcnt][maxEcnt].Tid = tagr->GetTid();
      npTag[minEcnt][maxEcnt].Status = 1;
    }  else {
      npTag[minEcnt][maxEcnt].Status = -1;
      return -1;
    }    
    npTag[minEcnt][maxEcnt].Sort = 1;
    return 0;
   } else {
     return -1 ;
   }
}


    
vector<Int_t> TNeutralPair::TagEnergyMatch( Bool_t (*testFunc) ( Double_t, vector<pointer>& ), 
					 vector<pointer>& pars, 
					 Int_t eBinMin, Int_t eBinMax, 
					 TH2D* hist ) {
  //
  // This function checks if there is a matching photon in 
  // tagger that satisfies the testFunc . 
  // The first parameter in testFunc is the photon beam energy, 
  // and the vector contains some parameter that testFunc will need 
  // to decide true or false (e.g. 4-vectors, cuts etc. )
  //
  vector<Int_t> matches;
  //  Int_t match = 0;
  if ( npReader != 0 && npOK ) {
    //
    // Loop to find the best photon comparing the photon RF time at the target center 
    // corrected for the neutrol z-position to the neutral vertex time
    //
    Int_t N_TAGR = npReader->GetNRows( "TAGR" );
    for( Int_t idxTAGR = 0; idxTAGR < N_TAGR; idxTAGR++ ) {
      TTAGRClass* tagr = static_cast<TTAGRClass*>( npReader->GetBankRow( "TAGR", idxTAGR ) );
      if( tagr != 0 && ( tagr->GetStat() == 15 || tagr->GetStat() == 7 ) && 
	  eBinMin < tagr->GetEid() && tagr->GetEid() < eBinMax && 
	  this->TidIsOK( tagr->GetTid() ) ) {
 	Double_t T_vert = npTagTimeOff + tagr->GetTagRF() + 
	  ( this->GetVtx().Z() - npTargetCenter ) / npClight; 
	Double_t Diff = TMath::Abs( npStartTime - T_vert ) ;
	if( testFunc( TNeutralPair::npMomCor.CorrectEgamma( tagr->GetEnergy() ), pars ) ) {
	  if( hist != 0 ) 
	    hist->Fill( tagr->GetTid(), npStartTime - T_vert );
	  if ( Diff  < 4*npTagrRes ) {
	    matches.push_back( idxTAGR + 1 ); 
	  }
	}	
      } 
    } 
  }  
  return matches;
}


void TNeutralPair::CorrectMomentum() {
  // 
  // Correct pi- and proton momentum for energy loss 
  // using eloss package. Then update the lambda momentum 
  // as well. 
  //

  // If Run Number changes, initialize ELoss package
  if ( npElossRunNo != npReader->GetRunNumber() ) {
    cout << "Switiching ELoss from run # " << npElossRunNo  << 
      " to run # " << npReader->GetRunNumber()  << endl;
    npElossRunNo = npReader->GetRunNumber() ;
    InitEloss( npElossRunNo );
  } 

  npNeg4Mom = TNeutralPair::ELossCorrect( npNeg4Mom, npVtx.Vect() );
  npPos4Mom = TNeutralPair::ELossCorrect( npPos4Mom, npVtx.Vect() );

  npNeg4Mom = TNeutralPair::npMomCor.CorrectMomentum( npNegId, npNeg4Mom );
  npPos4Mom = TNeutralPair::npMomCor.CorrectMomentum( npPosId, npPos4Mom );

  npMom = npNeg4Mom + npPos4Mom ;
  return;
}


TLorentzVector TNeutralPair::ELossCorrect( TLorentzVector lvP, TVector3 v3vtx ) {
  // Energy loss corrections 
  {
    Int_t iFlag = 2 ; // LD2 in the target
    Int_t iCell = 8 ; // EG3 cell
    
    vector4_t v4_p_meas = { lvP.T(), { lvP.X(), lvP.Y(), lvP.Z() } };
    vector3_t v3_vert = { v3vtx.X(), v3vtx.Y(), v3vtx.Z() };
    vector4_t v4_p_cor = c_momcor( v4_p_meas, lvP.M(), v3_vert, iFlag, iCell );
    
    return TLorentzVector( v4_p_cor.space.x, v4_p_cor.space.y, v4_p_cor.space.z, v4_p_cor.t );
  }
}



Int_t TNeutralPair::InitFit() {
  //
  // Initialize MINUIT by assigning title 
  // and setting verbose modes etc. 
  //
  if ( npOK ) {
    npFitGlob.KinFit.SetName( "NeutralFit" );
    npFitGlob.KinFit.SetTitle( "Neutral Fit" );

    //    npKinFit.Command( "SET STRategy 2" );
    npFitGlob.KinFit.Command( "SET NOWarnings"  );
    //    npKinFit.Command( "SET BAtch"  );

    npFitGlob.KinFit.SetPrintLevel( -1 );

  }
  return 0;
}



Int_t TNeutralPair::KinFit() {
  //
  // Perform a kinematic fit to lambda mass in two steps.
  // First find momenta for pion and proton and the opening angle 
  // using lambda mass constraint. Then using the obtained opening 
  // angle find the polar and azimuthal angles for proton and pion 
  // with angule fit.
  //
  if ( npOK ) {
    // Here we have to assign these global variable for the FCN functoins to access.
    // The only other alternative to the glabal variables is to pass 
    // numbers through unused parameters, which is messy too. 

    npFitGlob.MassPos = npMassPos;
    npFitGlob.MassNeg = npMassNeg;
    npFitGlob.MassNtr = npMassNtr;

    npFitGlob.Pos3Vec =  npPos4Mom.Vect() ;
    Double_t P1 = npFitGlob.Pos3Vec.Mag();
    Double_t The1 = npFitGlob.Pos3Vec.Theta();
    Double_t Phi1 = npFitGlob.Pos3Vec.Phi();
    Double_t dP1 = 0.015 * P1;
    Double_t dThe1 = 1.0 * TMath::DegToRad();
    Double_t dPhi1 = 1.0 * TMath::DegToRad();

    npFitGlob.Neg3Vec = npNeg4Mom.Vect() ;
    Double_t P2 = npFitGlob.Neg3Vec.Mag();
    Double_t The2 = npFitGlob.Neg3Vec.Theta();
    Double_t Phi2 = npFitGlob.Neg3Vec.Phi();
    Double_t dP2 = 0.015 * P2;
    Double_t dThe2 = 1.0 * TMath::DegToRad();
    Double_t dPhi2 = 1.0 * TMath::DegToRad();

    //    npFitGlob.Pos3Vec.Print();
    //    npFitGlob.Neg3Vec.Print();


    // First do fit with P1, P2 and opening angle The
    Double_t cst0 = cos( The1 ) * cos( The2 ) + sin( The1 ) * sin( The2 ) * cos( Phi2 - Phi1 );
    if ( TMath::Abs( cst0 ) > 1.0 ) cout << "Error a: Cosine greater than 1 : " << cst0 << endl;
    Double_t The = acos( cst0 );
    Double_t dThe = 1.0 * TMath::DegToRad();

    Int_t FitStat = FitMom( P1, P2, The, dP1, dP2, dThe );

     if ( FitStat == 0 ) {
      // Momentum Fit Converged
      Double_t P1_s, dP1_s, P2_s,  dP2_s, The_s, dThe_s ;
      
      npFitGlob.KinFit.GetParameter( 0, P1_s, dP1_s ); 
      npFitGlob.KinFit.GetParameter( 0, P2_s, dP2_s ); 
      npFitGlob.KinFit.GetParameter( 2, The_s, dThe_s ); 
      The_s = acos( cos( The_s ) );
      
      P2_s = P2fromTheP1( The_s, P1_s );
      if ( P2_s < 0 )  { 
	// Bad solution ?
	npFitMom.SetVectM( npMom.Vect(), npMassNtr ) ;	
	npXi2p = -1;
	npXi2 = -1;
	return -FitStat;
      } 
      
      Int_t nParDummy = 3;
      Double_t pars[] = { P1_s, P2_s, The_s };
      LbdFCN_P( nParDummy, static_cast<Double_t*>(0), npXi2p, static_cast<Double_t*>(pars), 1 ); 

      // Now fir the angles of proton and pi- to match the opening angle
      FitStat = FitAngles(  The1,  Phi1,  The2,  Phi2, 
			    dThe1, dPhi1, dThe2, dPhi2, The_s );

      if ( FitStat == 0 ) {
	// Angular Fit converget 
	Double_t  The1_s,  Phi1_s,  Phi2_s,  The2_s;
	Double_t dThe1_s, dPhi1_s, dPhi2_s, dThe2_s;
	npFitGlob.KinFit.GetParameter( 0, The1_s, dThe1_s );
	npFitGlob.KinFit.GetParameter( 1, Phi1_s, dPhi1_s );
	npFitGlob.KinFit.GetParameter( 2, The2_s, dThe2_s );
	npFitGlob.KinFit.GetParameter( 3, Phi2_s, dPhi2_s );

	The2_s = Th1fromThPh1Ph2( The1_s, Phi1_s, The2, Phi2_s, The_s );
	if ( The2_s < 0 ) {
	  //	  Bad solution ?
	  npFitMom.SetVectM( npMom.Vect(), npMassNtr ) ;	
	  npXi2a = -1;
	  npXi2 = -1;
	  return -FitStat;	  
	}

	Double_t pars[] = { The1_s, Phi1_s, The2_s, Phi2_s, The_s  };
	LbdFCN_Ang( nParDummy, static_cast<Double_t*>(0), npXi2a, static_cast<Double_t*>(pars), 1 ); 
	
	// Now set the npFitMom vector and add Xi2s
	TVector3 vProt;
	vProt.SetMagThetaPhi( P1_s, The1_s, Phi1_s );
	TVector3 vPion;
	vPion.SetMagThetaPhi( P2_s, The2_s, Phi2_s );
	
	npFitMom.SetVectM( ( vProt + vPion ), npMassNtr ) ; 
	npXi2 = npXi2p + npXi2a ;

      } else { 
	// Bad Migrad Result for angular fit
	npFitMom.SetVectM( npMom.Vect(), npMassNtr ) ;
	npXi2a = -1;
	npXi2 = -1;
	return FitStat;
      }
    } else {
      // Bad Migrad Result For momentum fit
      npFitMom.SetVectM( npMom.Vect(), npMassNtr ) ;
      npXi2p = -1;
      npXi2 = -1;
      return FitStat;
    }
    return FitStat;      
  } else 
    return -125;
}



Int_t TNeutralPair::FitMom( Double_t P1, Double_t P2, Double_t The, 
			    Double_t dP1, Double_t dP2,Double_t dThe ) {

  //
  // Do momentum fit with 2 free parameters
  //
  npFitGlob.KinFit.SetFCN( TNeutralPair::LbdFCN_P );
  
  npFitGlob.KinFit.DefineParameter( 0, "P1", P1, dP1, 0, 0 ) ;    
  npFitGlob.KinFit.DefineParameter( 1, "P2", P2, dP2, 0, 0 ) ;
  npFitGlob.KinFit.DefineParameter( 2, "The", The, dThe, 0, 0 ) ;
  npFitGlob.KinFit.DefineParameter( 3, "Dummy1", P2, dP2, 0, 0 ) ;
  npFitGlob.KinFit.DefineParameter( 4, "Dummy2", The, dThe, 0, 0 ) ;
  
  npFitGlob.KinFit.Release( 0 );
  npFitGlob.KinFit.FixParameter( 1 );
  npFitGlob.KinFit.Release( 2 ) ;
  npFitGlob.KinFit.FixParameter( 3 );
  npFitGlob.KinFit.FixParameter( 4 );
  
  return npFitGlob.KinFit.Migrad( );
}



Int_t TNeutralPair::FitAngles( Double_t The1,   Double_t Phi1,  Double_t The2,  Double_t Phi2,  
			       Double_t dThe1,Double_t  dPhi1, Double_t dThe2,  Double_t dPhi2, 
			       Double_t The0   ) {
  //
  // Do angular fit with 5 free parameters
  //
  npFitGlob.KinFit.SetFCN( TNeutralPair::LbdFCN_Ang );

  npFitGlob.KinFit.DefineParameter( 0, "The1", The1, dThe1, 0, 0 ) ;    
  npFitGlob.KinFit.DefineParameter( 1, "Phi1", Phi1, dPhi1, 0, 0 ) ;    
  npFitGlob.KinFit.DefineParameter( 2, "The2", The2, dThe2, 0, 0 ) ;    
  npFitGlob.KinFit.DefineParameter( 3, "Phi2", Phi2, dPhi2, 0, 0 ) ;    
  npFitGlob.KinFit.DefineParameter( 4, "The0", The0, 0.010, 0, 0 ) ;    

  npFitGlob.KinFit.Release( 0 );
  npFitGlob.KinFit.Release( 1 );
  npFitGlob.KinFit.FixParameter( 2 );
  npFitGlob.KinFit.Release( 3 );
  npFitGlob.KinFit.FixParameter( 4 );

  return npFitGlob.KinFit.Migrad();
}



void TNeutralPair::LbdFCN_P( Int_t& nPar, Double_t* Grad, Double_t& FunVal, Double_t* xVal, Int_t Flag ) {
  //
  // Function to minimize for momentum fit
  //
  Double_t P1 = npFitGlob.Pos3Vec.Mag();  
  Double_t The1 = npFitGlob.Pos3Vec.Theta();
  Double_t Phi1 = npFitGlob.Pos3Vec.Phi();
  Double_t dP1 = 0.015 * P1;
  
  Double_t P2 = npFitGlob.Neg3Vec.Mag(); 
  Double_t The2 = npFitGlob.Neg3Vec.Theta();
  Double_t Phi2 = npFitGlob.Neg3Vec.Phi();
  Double_t dP2 = 0.015 * P2;

  Double_t cst0 = cos( The1 ) * cos( The2 ) + sin( The1 ) * sin( The2 ) * cos( Phi2 - Phi1 );
  if ( TMath::Abs( cst0 ) > 1.0 ) cerr << "Error b: Cosine greater than 1 : " << cst0 << endl;
  Double_t The = acos( cst0 );
  Double_t dThe = 1.0 * TMath::DegToRad();

  Double_t P1_v = xVal[0];
  //  Double_t P2_v = xVal[1]; 
  Double_t The_v = xVal[2]; 

  // Bring The_v to [0, 180]
  The_v = acos( cos( The_v ) );
  
  Double_t P2_sol = P2fromTheP1( The_v, P1_v );
  if ( P2_sol < 0 ) {
    FunVal = +1.0e+34;
    return;
  }

  FunVal = 
    pow( P1_v   - P1, 2 )  / pow( dP1,  2 ) + 
    pow( P2_sol - P2, 2 )  / pow( dP2,  2 ) + 
    pow( The_v  - The, 2 ) / pow( dThe, 2) ;
  
  return;
}



void TNeutralPair::LbdFCN_Ang( Int_t& nPar, Double_t* Grad, Double_t& FunVal, Double_t* xVal, Int_t Flag ) {
  //
  // Function to minimize for angular fit.
  //
  Double_t The1 = npFitGlob.Pos3Vec.Theta();
  Double_t Phi1 = npFitGlob.Pos3Vec.Phi();
  Double_t dThe1 = 0.8 * TMath::DegToRad();
  Double_t dPhi1 = 1.3 * TMath::DegToRad();
  
  Double_t The2 = npFitGlob.Neg3Vec.Theta();
  Double_t Phi2 = npFitGlob.Neg3Vec.Phi();
  Double_t dThe2 = 0.8 * TMath::DegToRad();
  Double_t dPhi2 = 1.3 * TMath::DegToRad();

  Double_t The1_v = xVal[0];
  Double_t Phi1_v = xVal[1];
  //  Double_t The2_v = xVal[2];
  Double_t Phi2_v = xVal[3];
  Double_t The0   = xVal[4]; 
 

  // Bring phis to [0-360]
  if( Phi1_v > 0 ) 
    Phi1_v = fmod( Phi1_v, TMath::TwoPi() );
  if( Phi1_v < 0 ) 
    Phi1_v =  TMath::TwoPi() -  fmod( TMath::Abs( Phi1_v ), TMath::TwoPi() ) ;

  if( Phi2_v > 0 ) 
    Phi2_v = fmod( Phi2_v, TMath::TwoPi() );
  if( Phi2_v < 0 ) 
    Phi2_v =  TMath::TwoPi() -  fmod( TMath::Abs( Phi2_v ), TMath::TwoPi() ) ;

  // Bring The1 and Th2 to [0-180] 
  The1_v = acos( cos( The1_v ) );
  
  Double_t The2_s = Th1fromThPh1Ph2( The1_v, Phi1_v, The2, Phi2_v, The0 ); 
  if ( The2_s < 0 ) { 
    FunVal = +1.0e+34;
    return;
  }

  FunVal = 
    pow( The1_v  - The1, 2 ) / pow( dThe1, 2 ) + 
    pow( Phi1_v  - Phi1, 2 ) / pow( dPhi1, 2 ) + 
    pow( The2_s  - The2, 2 ) / pow( dThe2, 2 ) + 
    //    pow( xVal[2] - The2, 2 ) / pow( dThe2, 2 ) + 
    pow( Phi2_v  - Phi2, 2 ) / pow( dPhi2, 2 ) ; 

  return;

}



Double_t TNeutralPair::P2fromTheP1( Double_t The0, Double_t P1 ) {
  //
  // Function to calculate momentum of the pion from proton 
  // momentum and opening angle using the lambda mass constrainrt 
  // 
  Double_t M1 = npFitGlob.MassPos;
  Double_t M2 = npFitGlob.MassNeg;
  Double_t M0 = npFitGlob.MassNtr;  

  Double_t E1 = sqrt( P1*P1 + M1*M1 );
  Double_t cst = cos( The0 ) ;
  Double_t snt2 = 1.0 - cst*cst ;

  if ( ( M1*M1 + P1*P1*snt2 == 0 ) || 
       ( pow( M0*M0 - ( M1*M1 + M2*M2 ), 2 ) < ( 4*M2*M2 * ( M1*M1 + P1*P1*snt2 ) ) ) ) {
    return -1.0 ;
  }

  return ( 0.5 / ( M1*M1 + P1*P1*snt2 ) ) * 
    ( 
     ( ( M0*M0 - ( M1*M1 + M2*M2 ) ) * P1 * cst ) + 
     E1 * sqrt( pow( M0*M0 - ( M1*M1 + M2*M2 ), 2 ) - 4*M2*M2 * ( M1*M1 + P1*P1*snt2 )  )      
     );
}



Double_t TNeutralPair::Th1fromThPh1Ph2( Double_t The1, Double_t Phi1, 
					Double_t The2, Double_t Phi2  , Double_t The0 ) {
  //
  // Function to calculate pion polar angle from proton polar and azimuthal angles 
  // and pion azimuthal angle to match The0 opening angle. The The2 parameter 
  // is the original value so that we can pick the closest out of the 2 possible 
  // solutions for pion polar angle. This function can have 0, 1 or 2 solutions.
  //
  Double_t cst0 = cos( The0 );
  Double_t cst1 = cos( The1 );
  Double_t cst2 = cos( The2 );
  Double_t snt1 = sin( The1 );
  Double_t cos_dphi = cos( Phi2 - Phi1 ) ;
  
  if ( cst1 == 0.0 && cos_dphi == 0.0 ) {
    // Here any theta is good, so just use the original one.
    // Should return to avoid division by zero
    return The2 ;
  } else { 
    if ( ( -cst0*cst0 + cst1*cst1 + snt1*snt1 * cos_dphi*cos_dphi ) < 0 ) {
      // This is when the plane defined by phi2 is two far from 
      // the (the1, phi1) direction and cannot be reached by The0 cone
      //     cout << "No solution found " << endl;
      return -2;
    }
    // There are two solutions (sometimes just one ) 
    // pick the one closest to the original
    Double_t sol_plus  = ( cst0 * cst1  + snt1  * cos_dphi * 
		       sqrt( -cst0*cst0 + cst1*cst1 + snt1*snt1 * cos_dphi*cos_dphi ) ) / 
      ( cst1*cst1 + snt1*snt1 * cos_dphi*cos_dphi ) ;

    Double_t sol_minus = ( cst0 * cst1  - snt1  * cos_dphi * 
		       sqrt( -cst0*cst0 + cst1*cst1 + snt1*snt1 * cos_dphi*cos_dphi ) ) / 
      ( cst1*cst1 + snt1*snt1 * cos_dphi*cos_dphi ) ;
    
    if ( TMath::Abs ( sol_plus ) <= 1 && TMath::Abs( sol_minus ) <= 1.0 ) {
      if (TMath::Abs( sol_plus - cst2 ) > TMath::Abs( sol_minus - cst2 ) ) 
	return acos( sol_minus ) ;
      else 
	return acos( sol_plus );
    }
    else { 
      if( TMath::Abs ( sol_minus ) <= 1 ) {
	return acos( sol_minus ) ;
      } 
      else { 
	if ( TMath::Abs ( sol_plus ) <= 1 ) {
	  return acos( sol_plus );
	} 
	else {
	  cout << "Bad values for cosines " << sol_plus << " " << sol_minus << endl; 
	  return -1;
	}
      }
    }
  }
}
