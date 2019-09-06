#ifndef _TTBIDClass_
#define _TTBIDClass_

#include <iostream>
using namespace std;

#include "TObject.h"
#include "TString.h"

class TTBIDClass: public TObject{

public:
	Int_t track;		// Track (index to TBTR) or zero if neutral
	Int_t sec;		// Sector track is in
	Float_t beta;		// Beta of the track in units of c
	Float_t vtime;		// vertex time of track

	Int_t sc_stat;		// status of hit matching to SC: see sc.h
	Int_t sc_id;		// Pointer to SCRC bank
	Float_t sc_time;	// SC calibrated time for this track (ns)
	Float_t sc_qual;	// quality of match for SC
	Float_t sc_vtime;	// time at vertex for SC(ns)
	Float_t sc_beta;	// Beta calculated from TOF from SC

	Int_t cc_stat;		// status of hit matching to CC: see sc.h
	Int_t cc_id;		// pointer to CC01 bank
	Float_t cc_time;	// CC calibrated time for this track (ns)
	Float_t cc_qual;	// quality of match for CC
	Float_t	cc_vtime;	// time at vertex for CC(ns)
	Float_t cc_beta;	// Beta as calculated by the EC

	Int_t ec_stat;		// status of hit matching to ec: see sc.h
	Int_t ec_id;		// Pointer to ECHB bank
	Float_t ec_time;	// EC calibrated time for this track (ns)
	Float_t ec_qual;	// EC quality factor
	Float_t ec_vtime;	// time at vertex for EC(ns)
	Float_t ec_beta;	// Beta as calculated by the EC

	Int_t st_stat;		// status of hit matching to ST
	Int_t st_id;		// Pointer to STR bank
	Float_t st_time;	// ST calibrated time for this track (ns)
	Float_t st_qual;	// ST quality factor
	Float_t st_vtime;	// time at vertex for ST(ns)
	Float_t st_beta;	// Beta as calculated by the ST

	Int_t lac_stat;		// status of hit matching to LAC
	Int_t lac_id;		// Pointer to EC1R bank
	Float_t lac_time;	// LAC calibrated time for this track (ns)
	Float_t lac_qual;	// LAC quality factor
	Float_t lac_vtime;	// time at vertex for LAC(ns)
	Float_t lac_beta;	// Beta as calculated by the LAC

public:
	TTBIDClass(){};
	~TTBIDClass(){};
	TTBIDClass(TTBIDClass *TmpTBID); // copy constructor.

	inline Int_t GetTrack() { return track; }
	inline Int_t GetSec()   { return sec; }
	inline Float_t GetBeta() { return beta; }
	inline Float_t GetVtime() { return vtime; }

	inline Int_t GetSc_stat() { return sc_stat; }
	inline Int_t GetSc_id()   { return sc_id; }
	inline Float_t GetSc_time() { return sc_time; }
	inline Float_t GetSc_qual() { return sc_qual; }
	inline Float_t GetSc_vtime() { return sc_vtime; }
	inline Float_t GetSc_beta()  { return sc_beta; }

	inline Int_t GetCc_stat()    { return cc_stat; }
	inline Int_t GetCc_id()      { return cc_id; }
	inline Float_t GetCc_time()  { return cc_time; }
	inline Float_t GetCc_qual()  { return cc_qual; }
	inline Float_t GetCc_vtime() { return cc_vtime; }
	inline Float_t GetCc_beta()  { return cc_beta; }

	inline Int_t GetEc_stat()    { return ec_stat; }
	inline Int_t GetEc_id()      { return ec_id; }
	inline Float_t GetEc_time()  { return ec_time; }
	inline Float_t GetEc_qual()  { return ec_qual; }
	inline Float_t GetEc_vtime() { return ec_vtime; }
	inline Float_t GetEc_beta()  { return ec_beta; }

	inline Int_t GetSt_stat()    { return st_stat; }
	inline Int_t GetSt_id()      { return st_id; }
	inline Float_t GetSt_time()  { return st_time; }
	inline Float_t GetSt_qual()  { return st_qual; }
	inline Float_t GetSt_vtime() { return st_vtime; }
	inline Float_t GetSt_beta()  { return st_beta; }

	inline Int_t GetLac_stat()   { return lac_stat; }
	inline Int_t GetLac_id()     { return lac_id; }
	inline Float_t GetLac_time() { return lac_time; }
	inline Float_t GetLac_qual() { return lac_qual; }
	inline Float_t GetLac_vtime() { return lac_vtime; }
	inline Float_t GetLac_beta()  { return lac_beta; }

	void Print(void);   // Print info on class.

	ClassDef(TTBIDClass,1) // Class for accessing the TBER bank: Tracking error information.
};
#endif
