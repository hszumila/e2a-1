
#include "TItemValue.h"

ClassImp(TItemValue)


TItemValue::TItemValue(){
// default constructor
//
}

TItemValue::~TItemValue(){
// default destructure
//
}

TItemValue::TItemValue(Int_t fpItemValueId, Int_t    fpItemId,Int_t    fpMinRun, Int_t    fpMaxRun, const char *fpAuthor, const char *fpTime, const char *fpComment){
// Constructor Sets Variables 
//
fItemValueId =  fpItemValueId;
 fItemId = fpItemId;
  fMinRun = fpMinRun;
   fMaxRun = fpMaxRun;
    fAuthor = fpAuthor; 
     fTime = fpTime;
      fComment = fpComment;

}

void TItemValue::Set(Int_t fpItemValueId, Int_t    fpItemId,Int_t    fpMinRun, Int_t    fpMaxRun, const char *fpAuthor, const char *fpTime, const char *fpComment){
// Sets Variables 
//

fItemValueId = fpItemValueId;
 fItemId = fpItemId;
  fMinRun = fpMinRun;
   fMaxRun = fpMaxRun;
    fAuthor = fpAuthor; 
     fTime = fpTime;
      fComment = fpComment;

}

