#ifndef _TTBTRClass_
#define _TTBTRClass_

#include <iostream>
using namespace std;

#include "TObject.h"
#include "TString.h"
#include "TMatrix.h"

class TTBTRClass: public TObject{

public:
	Float_t x;		// x
	Float_t y;		// y 'Vertex' position {x,y,z}
	Float_t z;		// z
	Float_t px;		// Px
	Float_t py;		// Py  momentum at 'vertex' {Px,Py,Pz}
	Float_t pz;		// Pz
	Float_t q;		// charge   (straight tracks: set q=0.)
	Float_t chi2;		// Chisquare for this Track
	Int_t itr_sec;		// Trackno_in_Sector + Sector*100
	Int_t itr_hbt;		// Trackno in HBTR for this track

public:
	TTBTRClass(){};
	~TTBTRClass(){};
	TTBTRClass(TTBTRClass *TmpTBTR); // copy constructor.

	inline Float_t GetX() { return x; }
	inline Float_t GetY() { return y; }
	inline Float_t GetZ() { return z; }
	inline Float_t GetPx()  { return px; }
	inline Float_t GetPy()  { return py; }
	inline Float_t GetPz()  { return pz; }
	inline Float_t GetCharge() { return q; }
	inline Float_t GetChi2()  { return chi2; }
	inline Int_t GetItr_sec()  { return itr_sec; }
	inline Int_t GetItr_hbt()  { return itr_hbt; }

	void Print(void);   // Print info on class.

	ClassDef(TTBTRClass,1) // Class for accessing the TBER bank: Tracking error information.
};
#endif
