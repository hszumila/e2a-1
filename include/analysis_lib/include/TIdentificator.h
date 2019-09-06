#ifndef TIDENTIFICATOR_H
#define TIDENTIFICATOR_H

#include "TClasTool.h"


class TIdentificator {
public:
    explicit TIdentificator(TClasTool* CT = 0);
    ~TIdentificator();

    // HEADER bank
    Float_t NEvent();                             // inline
    Float_t FC();                                 // inline
    Float_t FCG();                                // inline - ADDED BY Barak - 05.12.2015 
    Float_t STT();                                // inline
    Float_t NRun();                               // inline - ADDED By Barak - 11.20.2017

    // EVNT bank
    Double_t Betta(Int_t k);                      // inline
    Double_t Id(Int_t k, Bool_t kind = 0);        // inline
    Double_t Charge(Int_t k);                     // inline
    Double_t Px(Int_t k, Bool_t kind = 0);        // inline
    Double_t Py(Int_t k, Bool_t kind = 0);        // inline
    Double_t Pz(Int_t k, Bool_t kind = 0);        // inline
    Double_t X(Int_t k);                          // inline
    Double_t Y(Int_t k);                          // inline
    Double_t Z(Int_t k, Bool_t kind = 0);         // inline
    Double_t Z_Corrected(Int_t k, TString particle);
    Int_t StatSC(Int_t k);                        // inline
    Int_t StatCC(Int_t k);                        // inline
    Int_t StatDC(Int_t k);                        // inline
    Int_t StatEC(Int_t k);                        // inline
    Double_t Status(Int_t k);                     // inline
    Int_t CCPart();                               // inline - ADDED BY Barak - 05.21.2015 
    Int_t DCPart();                               // inline - ADDED BY Barak - 05.21.2015 
    Int_t ECPart();                               // inline - ADDED BY Barak - 05.21.2015 
    Int_t SCPart();                               // inline - ADDED BY Barak - 05.21.2015
    Int_t gPart();                                // inline - ADDED BY Barak - 05.21.2015

    // CCPB
    Double_t Nphe(Int_t k);                       // inline
    Double_t Chi2CC(Int_t k);                     // inline
    Double_t CCStatus(Int_t k);                   // inline
    Double_t PathCC(Int_t k);                     // inline - ADDED By Barak - 12.06.2017
    Double_t TimeCC(Int_t k);                     // inline - ADDED By Barak - 12.06.2017
    
    // DCPB
    Double_t DCStatus(Int_t k);                   // inline

    // ECPB
    Double_t Etot(Int_t k);                       // inline
    Double_t Ein(Int_t k);                        // inline
    Double_t Eout(Int_t k);                       // inline
    Double_t EC_path(Int_t k);                    // inline
    Double_t EC_time(Int_t k);                    // inline
    Double_t EC_U(Int_t k);                       // inline
    Double_t EC_V(Int_t k);                       // inline
    Double_t EC_W(Int_t k);                       // inline
    Double_t EC_X(Int_t k);                       // inline
    Double_t EC_Y(Int_t k);                       // inline
    Double_t EC_Z(Int_t k);                       // inline

    Double_t ECStatus(Int_t k);                   // inline

    // SCPB
    Double_t PathSC(Int_t k);                     // inline
    Double_t TimeSC(Int_t k);                     // inline
    Double_t EdepSC(Int_t k);                     // inline
    Double_t SCStatus(Int_t k);                   // inline
    Double_t SCpdht(Int_t k);			  // inline - ADDED BY MISHA - 07.01.2014

    // Derived observables
    Double_t Momentum(Int_t k, Bool_t kind = 0);
    Double_t Moment(Int_t k, Bool_t kind = 0);    // Deprecated
    Double_t Mass2(Int_t k);
    Double_t ThetaLab(Int_t k, Bool_t kind = 0);
    Double_t PhiLab(Int_t k, Bool_t kind = 0);
    Double_t ThetaVirtLab(Bool_t kind = 0);
    Double_t PhiVirtLab(Bool_t kind = 0);
    Double_t ThetaPQ(Int_t k, Bool_t kind = 0);   // See Barak's Comments
    Double_t PhiPQ(Int_t k, Bool_t kind = 0);
    Double_t CosThetaPQ(Int_t k, Bool_t kind = 0);
    Double_t PTrans2PQ(Int_t k, Bool_t kind = 0);
    Double_t PLong2PQ(Int_t k, Bool_t kind = 0);
    Int_t Sector(Int_t k, Bool_t kind = 0);

    // Kinematic variables
    Double_t Q2(Bool_t kind = 0);
    Double_t W(Bool_t kind = 0);
    Double_t Nu(Bool_t kind = 0);
    Double_t Xb(Bool_t kind = 0);
    Double_t Yb(Bool_t kind = 0);
    Double_t ZhPi(Int_t k, Double_t Mass, Bool_t kind = 0);


    // Correction functions
    Double_t TimeCorr4(Double_t mass, Int_t k);

    // Particle Identification cuts
    TString GetCategorization(Int_t k);
    TString* GetCategorization();
    void PrintCategorization();
    void PrintCategorization(TString* partIds);
    bool ProtonPID(Double_t time, Double_t p);
    Double_t electronPID_b(Double_t p);
    //bool ProtonPID_b(Double_t edep, Double_t p);  // ADDED BY Barak - 12.07.2015


    // Fiducial Cut
    // Pi-Minus Cuts Added by Barak - 10.23.2017
    Double_t FidTheta(Int_t k, Bool_t kind = 0);
    Double_t FidThetaMin();
    Double_t FidThetaMinPiPlus(Int_t k);
    Double_t FidThetaMinPiMinus(Int_t k);
    Double_t FidFunc(Int_t side, Int_t param);
    Double_t FidFuncPiPlus(Int_t side, Int_t param, Int_t k);
    Double_t FidFuncPiMinus(Int_t side, Int_t param, Int_t k);
    Double_t FidPhi(Int_t k, Bool_t kind = 0);
    Double_t FidPhiMin();
    Double_t FidPhiMax();
    Double_t FidPhiMinPiPlus(Int_t k);
    Double_t FidPhiMaxPiPlus(Int_t k);
    Double_t FidPhiMinPiMinus(Int_t k);
    Double_t FidPhiMaxPiMinus(Int_t k);
    Bool_t FidCheckCut();
    Bool_t FidCheckCutPiPlus(Int_t k);
    Bool_t FidCheckCutPiMinus(Int_t k);
    Int_t FidSector(Int_t k, Bool_t kind = 0);


    //Target methods.
    Int_t ElecVertTarg();
    //Int_t ElecVertTargCorrected();
    Bool_t PionVertTarg(Int_t k);


    //New Corrections - ADDED BY Barak - 08.20.2015, 05.06.2016
    Double_t Q2_Corrected(Bool_t kind = 0);
    Double_t Nu_Corrected(Bool_t kind = 0);
    Double_t W_Corrected(Bool_t kind = 0);
    Double_t Xb_Corrected(Bool_t kind = 0);
    Double_t Yb_Corrected(Bool_t kind = 0);
    Double_t CosThetaPQ_Corrected(Int_t k, Bool_t kind = 0);

    Double_t Momentum_Corrected(Int_t k);
    Double_t Theta_Corrected(Int_t k); //Angle in degrees
    //Double_t Z_Good(Int_t k);
    //Double_t Z_Better(Int_t k);

private:
    const Double_t kEbeam;    // The energy of incoming electron beam
    const Double_t kEbeam_arc;// The energy of incoming electron beam (from Hall A Arc) - ADDED BY Barak - 08.20.2015 
    const Double_t kMpi;      // The mass of the pion
    const Double_t kGOOD;     // The key for the exceptions (should be improved to avoid it at all !!!)

    TClasTool *fCT;           // Pointer to the main ClasTool object
    
    TEVNTClass *fEVNT;        // Pointer to the EVNT object
    TGSIMClass *fGSIM;        // Pointer to the GSIM object
    TCCPBClass *fCCPB;        // Pointer to the CCPB object
    TECPBClass *fECPB;        // Pointer to the ECPB object
    TSCPBClass *fSCPB;        // Pointer to the SCPB object
    TDCPBClass *fDCPB;        // Pointer to the DCPB object

    TString* fPartIds;        // Array with the categories of the particles belonging to an event.
};

#include "TIdentificator.icc"

#endif
