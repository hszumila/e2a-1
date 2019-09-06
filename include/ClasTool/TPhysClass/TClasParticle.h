// Autogenerated Class (Header File)
// Author : G.Gavalian
// Date   : Tue Jul 21 08:39:39 EDT 2009
//

#ifndef __TClasParticle__
#define __TClasParticle__
#include <iostream>
#include <TROOT.h>
#include <TVector3.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <sstream>
#include <iomanip>

using namespace std;

class TClasParticle : public TObject {

private:

  int     pdgID;
  double  pdgMass;
  double  tofMass;
  int     partCharge;

  TLorentzVector  fPartVector;
  TVector3        fPartVertex;
  Int_t           fPartStatus;


  string  vectorToString(TVector3 v3, const char *_pre = "(",
			 const char *_post = ")");
  string  particleToStringLund();

public:

TClasParticle ();
~TClasParticle ();
 TClasParticle(const TClasParticle &part);


 void   Reset();
 void   SetId(int id);

 void   SetStatus(int stat){
   if(stat>0) fPartStatus = 1; else fPartStatus = 0;
 };

 void   SetTofMass(double mass){tofMass = mass;};
 void   SetPdgMass(double mass){pdgMass = mass;};
 void   SetVector(TLorentzVector vect){fPartVector = vect;};
 void   SetVertex(TVector3  vect){fPartVertex = vect;};
 void   SetPartById(int id, TVector3 p3, TVector3 v3);
 
 TLorentzVector   GetVector(){return fPartVector;};
 TVector3         GetVertex(){ return fPartVertex;};
 Double_t         GetPdgMass(){return pdgMass;};
 Double_t         GetTofMass(){return tofMass;};
 Int_t            GetPdgId(){return pdgID;};
 Int_t            GetStatus(){return fPartStatus;};
 TString          ToString();
 TClasParticle &operator=(const TClasParticle &part);



ClassDef(TClasParticle,1)


};
#endif
