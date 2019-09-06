////////////////////////////////////////////////////////////////////////
// File TNeutralPair.h
//
// Initial Author: Hovanes Egiyan  UNH  11/21/05
//
// Clas to find best Neutral pair 
// 
////////////////////////////////////////////////////////////////////////
 
#ifndef _TNeutralPair_
#define _TNeutralPair_
#include <iostream>
#include <map>
#include <cmath>
#include <vector>
using namespace std;
#include "TVector3.h" 
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TH2.h"
    
#include "THEADERClass.h"
#include "TEVNTClass.h"
#include "TTAGRClass.h"
#include "TDCPBClass.h"
#include "TSCPBClass.h"
#include "TSTPBClass.h"
#include "TECPBClass.h"
#include "TCCPBClass.h"
#include "TDSTReader.h"
#include "TClasTool.h"
#include "TMomCorEG3.h"

     
extern "C" {
#include "kinematics.h"
#include "eloss.h" 
} 

#define N_ECNT (767)

typedef void* pointer;      // Had to do it because of stupid rootcint 



 
class TNeutralPair {

  // private: 
 public: 

  class FitGlobType {  // A structure to be used by MINUIT
  public :
    TMinuit  KinFit;      // Global variable for MINUIT fitting
    Bool_t   MinuitFlag;  // Variable to indicate if MINUIT was inited
    TVector3 Pos3Vec;  // Global var for Pos 3-vec
    TVector3 Neg3Vec;  // Global var for Neg 3-vec
    Double_t MassPos;  // Gloal var for Positive Mass
    Double_t MassNeg;  // Global var for Negative Mass
    Double_t MassNtr; // Global var fot Neutral Mass
  } ;
  
  class TTagInfo {
  public: 
    Int_t Status;
    Int_t Sort;
    Int_t Eid;
    Int_t Tid;
    Double_t Energy;
    Double_t Time;

    TTagInfo( Int_t status = 0, Int_t sort = 0, 
	      Int_t eid    = 0, Int_t tid  = 0, 
	      Double_t energy = 0, Double_t time = 0 ) 
      : Status(status), Sort(sort), Eid(eid), 
      Tid(tid), Energy(energy), Time(time ) {}

    TTagInfo( const TNeutralPair::TTagInfo& tag ) 
      : Status(tag.Status), Sort(tag.Sort), Eid(tag.Eid), 
      Tid(tag.Tid), Energy(tag.Energy), Time(tag.Time ) {}
    
    TTagInfo operator=( TNeutralPair::TTagInfo tag ) ;
  };

  Bool_t npOK;

  Double_t npStartTime;   // start time from the neutral 

  Bool_t npKinFitFlag;    // Flag to show if kinematic flag was requested
  Bool_t npElossFlag;     // Flag to indicate if Mom. corrections done

  Int_t npPosId;         // ID for positive particle
  Int_t npNegId;         // ID for negative particle 
  Int_t npNtrId;         // ID for the combined particle

  Int_t npBestPosIdx;    // EVNT index for positive 
  Int_t npBestNegIdx;    // EVNT index for negative
 
  Double_t npBestMass;   // Mass for the best pair
  Double_t npBestDOCA;   // DOCA for the best pair
  Double_t npBestVtxTimeDiff;  // Time difference for the best pair

  Double_t npXi2;         // Xi2 from the overall fit
  Double_t npXi2p;        // Xi2 part from momentum fit
  Double_t npXi2a;        // Xi2 from angular fit

  TDSTReader* npReader;   // reader to use

  TLorentzVector npNeg4Mom; // 4-mom of negative
  TLorentzVector npPos4Mom; // 4-mom of positve

  TLorentzVector  npMom;  // 4-mom of the neutral
  TLorentzVector  npVtx;  // 4-vertex of the neutral  
  
  TLorentzVector  npFitMom; // momentum after fit

  TVectorD* npBeta;      // Vector for betas for all tracks in EVNT

  Double_t npMassPos;   // Mass of the positive
  Double_t npMassNeg;   // Mass of the negative 
  Double_t npMassNtr;   // Mass of combined system 

  Double_t npNeutralPeak; // Position of the neutral peak w/o corections
  map<Int_t, map<Int_t, TTagInfo> >  npTag;

  vector<bool> npTidIsOK;   // Vector for T-d statuses 

  static FitGlobType npFitGlob;  // Must be global to be used by MINUIT
  static TDatabasePDG npDbPDG;

  static Double_t npClight; // Speed of light 

  static Double_t npSigmaMass; // Mass resolution
  static Double_t npSigmaDOCA; // DOCA resolution
  static Double_t npSigmaTime; // Time resolution
 
  static Double_t npMassCut;   // scale factor for mass cut 
  static Double_t npDOCACut;   // scale factor for DOCA cut
  static Double_t npTimeCut;   // scale factor for Time cut

  static Double_t npTagrRes;  // Tagger resolition
  static Double_t npTargetCenter; // Target center 
  static Double_t npTagrErgShift; // Tagger energy shift

  static Double_t npTagTimeOff;   // Time offset between pair time and RF time
  static Double_t npSimLbdTimeOff; // Time offset between negative and positive for GSIM 

  static Double_t npMinMomPos;    // Minimum momentum for a positive tack to have
  static Double_t npMinMomNeg;    // Minimum momentum for a negative tack to have

  static Int_t    npElossRunNo;  // Run number for the last eloss call




  Int_t     FindBestNeutral();
  Double_t  FindFlightTime();
  void      FindBetas();
  Double_t  FindStartTime();
  
  void      CorrectMomentum();

  Int_t     InitFit();
  Int_t     KinFit();

  Int_t FitMom( Double_t P1, Double_t P2, Double_t The, 
		Double_t dP1, Double_t dP2,Double_t dThe );
  Int_t FitAngles( Double_t The1,   Double_t Phi1,  Double_t The2,  Double_t Phi2,  
			      Double_t dThe1,Double_t  dPhi1, Double_t dThe2,  Double_t dPhi2, 
		   Double_t The0   );
  
  static Double_t P2fromTheP1( Double_t The0, Double_t P1 );
  static Double_t Th1fromThPh1Ph2( Double_t The1, Double_t Phi1, 
						Double_t The2, Double_t Phi2  , Double_t The0 );


 public:
  static TMomCorEG3  npMomCor;

  TNeutralPair( TDSTReader* reader,  Bool_t eLossFlag = 0, Bool_t fitFlag = 0, 
		Int_t ntrID = 3122, Int_t posId = 2212, Int_t negId = -211, 
		Double_t mPeak = -10000 );
  ~TNeutralPair( );
  TNeutralPair( TNeutralPair& pair );

  TNeutralPair& operator=( TNeutralPair& pair );

  vector<Int_t> TagEnergyMatch( Bool_t (*testFunc) ( Double_t, vector<pointer>& ), 
			     vector<pointer>& pars, 
			     Int_t eBinMin = 1, Int_t eBinMax = N_ECNT, 
			     TH2D* hist = 0 );

  Double_t GetMatchTime( TVector3 p, TVector3 v );

  static TLorentzVector ELossCorrect( TLorentzVector lvP, TVector3 v3vtx );
  
  Int_t    TagTimeSort( Int_t minEcnt = 1, Int_t maxEcnt = N_ECNT  );
  Double_t GetTagTime( Int_t minEcnt = 1, Int_t maxEcnt = N_ECNT ) ;
  Double_t GetTagEnergy( Int_t minEcnt = 1, Int_t maxEcnt = N_ECNT ) ;
  Int_t    GetErgBin( Int_t minEcnt = 1, Int_t maxEcnt = N_ECNT );
  Int_t    GetTimeBin( Int_t minEcnt = 1, Int_t maxEcnt = N_ECNT );

  static void  LbdFCN_P( Int_t& nPar, Double_t* Grad, Double_t& FunVal, Double_t* xVal, Int_t Flag ) ;
  static void  LbdFCN_Ang( Int_t& nPar, Double_t* Grad, Double_t& FunVal, Double_t* xVal, Int_t Flag ) ;

  static TVector3 FindNeutralVtx( TVector3 v1, TVector3 p1, 
			  TVector3 v2, TVector3 p2, Double_t *dist );

  inline Bool_t IsOK() 
    { return npOK ; }

  inline Bool_t GetElossFlag() 
    { return npElossFlag; } 
  
  inline Bool_t GetKinFitFlag() 
    { return npKinFitFlag; }


  inline Double_t GetStartTime() 
    { return npStartTime; }

  inline Double_t GetBestMass() 
    { return npBestMass ; }
  inline Double_t GetBestDOCA() 
    { return npBestDOCA ; }
  inline Double_t GetBestTimeDiff() 
    { return npBestVtxTimeDiff ; }


  inline Int_t GetPosIdx() 
    { return npBestPosIdx ; }
  inline Int_t GetNegIdx() 
    { return npBestNegIdx ; }

  inline TLorentzVector GetPos4Mom() 
    { return npPos4Mom; }
  inline TLorentzVector GetNeg4Mom()
    { return npNeg4Mom; }
  
  inline TLorentzVector GetMom() 
    { return npMom; }
  inline TLorentzVector GetVtx() 
    { return npVtx; }

  inline TLorentzVector GetFitMom()
    { return npFitMom; }

  inline TDSTReader* GetReaderPtr() 
    { return npReader; }

  inline Double_t GetXi2() 
    { return npXi2 ; }

  inline Double_t GetMomXi2() 
    { return npXi2p ; }

  inline Double_t GetAngXi2() 
    { return npXi2a ; }

  inline Double_t GetBeta( Int_t iTrk ) 
    { if ( npBeta != 0 && iTrk < npBeta->GetNoElements() )
	return (*npBeta)[iTrk] ;
      else 
	return 0.0 ; 
    }

  inline Double_t SetNeutralPeak( Double_t peak ) 
    { return ( npNeutralPeak = peak ); }
  inline  Double_t GetNeutralPeak( ) 
    { return npNeutralPeak ; }


  inline void SetTidStatus( vector<bool>& tidIsOK )
  { npTidIsOK = tidIsOK; return; }

  inline Bool_t TidIsOK ( unsigned int tID ) {
    if( ( tID <= npTidIsOK.size() ) && npTidIsOK[tID-1] ) 
      return true;
    else 
      return false;
  }

  inline static Double_t SetTagrRes( Double_t res ) 
  { return ( npTagrRes = res ) ; } 
  inline static Double_t GetTagrRes( Double_t res ) 
  { return  npTagrRes ; } 
  

  inline static Double_t SetSigmaMass( Double_t sigma ) 
    { return ( TNeutralPair::npSigmaMass = sigma ); }
  inline static Double_t GetSigmaMass( ) 
    { return TNeutralPair::npSigmaMass ; }

  inline static Double_t SetSigmaDOCA( Double_t sigma ) 
    { return ( TNeutralPair::npSigmaDOCA = sigma ); }
  inline static Double_t GetSigmaDOCA( ) 
    { return TNeutralPair::npSigmaDOCA ; }
  

  inline static Double_t SetSigmaTime( Double_t sigma ) 
    { return ( TNeutralPair::npSigmaTime = sigma ); }
  inline static Double_t GetSigmaTime( )  
    { return TNeutralPair::npSigmaTime ; }

  
  inline static Double_t SetMassCut( Double_t cut ) 
    { return ( TNeutralPair::npMassCut = cut ); }
  inline static Double_t GetMassCut( ) 
    { return TNeutralPair::npMassCut ; }


  inline static Double_t SetDOCACut( Double_t cut ) 
    { return ( TNeutralPair::npDOCACut = cut ); }
  inline static Double_t GetDOCACut( ) 
    { return TNeutralPair::npDOCACut ; }

  inline static Double_t SetTimeCut( Double_t cut ) 
    { return ( TNeutralPair::npTimeCut = cut ); }
  inline static Double_t GetTimeCut( ) 
    { return TNeutralPair::npTimeCut ; }

  static inline Int_t GetElossRunNo() 
    { return npElossRunNo ; }

  static inline Double_t SetMinMomPos( Double_t momVal )
    { return ( TNeutralPair::npMinMomPos = momVal ); }
  static inline Double_t GetMinMomPos()
    { return TNeutralPair::npMinMomPos ; }

  static inline Double_t SetMinMomNeg( Double_t momVal )
    { return ( TNeutralPair::npMinMomNeg = momVal ) ; }
  static inline Double_t GetMinMomNeg()
    { return TNeutralPair::npMinMomNeg ; }

  static inline Double_t SetTagRes( Double_t resol ) 
    { return ( TNeutralPair::npTagrRes = resol ) ; }
  static inline Double_t GetTagRes() 
    { return TNeutralPair::npTagrRes ; }

};




#endif
