//
// File generated by /u/apps/root/5.34.13/root/bin/rootcint at Fri Jan 20 19:15:40 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL6dITSCPBClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TSCPBClassDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TSCPBClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TSCPBClass(void *p = 0);
   static void *newArray_TSCPBClass(Long_t size, void *p);
   static void delete_TSCPBClass(void *p);
   static void deleteArray_TSCPBClass(void *p);
   static void destruct_TSCPBClass(void *p);
   static void streamer_TSCPBClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TSCPBClass*)
   {
      ::TSCPBClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TSCPBClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TSCPBClass", ::TSCPBClass::Class_Version(), "./TSCPBClass.h", 24,
                  typeid(::TSCPBClass), DefineBehavior(ptr, ptr),
                  &::TSCPBClass::Dictionary, isa_proxy, 0,
                  sizeof(::TSCPBClass) );
      instance.SetNew(&new_TSCPBClass);
      instance.SetNewArray(&newArray_TSCPBClass);
      instance.SetDelete(&delete_TSCPBClass);
      instance.SetDeleteArray(&deleteArray_TSCPBClass);
      instance.SetDestructor(&destruct_TSCPBClass);
      instance.SetStreamerFunc(&streamer_TSCPBClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TSCPBClass*)
   {
      return GenerateInitInstanceLocal((::TSCPBClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TSCPBClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TSCPBClass::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TSCPBClass::Class_Name()
{
   return "TSCPBClass";
}

//______________________________________________________________________________
const char *TSCPBClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSCPBClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TSCPBClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSCPBClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TSCPBClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSCPBClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TSCPBClass::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSCPBClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TSCPBClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class TSCPBClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Scpdht;
      R__b >> Edep;
      R__b >> Time;
      R__b >> Path;
      R__b >> Chi2sc;
      R__b >> Status;
      R__b.CheckByteCount(R__s, R__c, TSCPBClass::IsA());
   } else {
      R__c = R__b.WriteVersion(TSCPBClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Scpdht;
      R__b << Edep;
      R__b << Time;
      R__b << Path;
      R__b << Chi2sc;
      R__b << Status;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TSCPBClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TSCPBClass.
      TClass *R__cl = ::TSCPBClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Scpdht", &Scpdht);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Edep", &Edep);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Time", &Time);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Path", &Path);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Chi2sc", &Chi2sc);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Status", &Status);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TSCPBClass(void *p) {
      return  p ? new(p) ::TSCPBClass : new ::TSCPBClass;
   }
   static void *newArray_TSCPBClass(Long_t nElements, void *p) {
      return p ? new(p) ::TSCPBClass[nElements] : new ::TSCPBClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_TSCPBClass(void *p) {
      delete ((::TSCPBClass*)p);
   }
   static void deleteArray_TSCPBClass(void *p) {
      delete [] ((::TSCPBClass*)p);
   }
   static void destruct_TSCPBClass(void *p) {
      typedef ::TSCPBClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TSCPBClass(TBuffer &buf, void *obj) {
      ((::TSCPBClass*)obj)->::TSCPBClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::TSCPBClass

/********************************************************
* dicts/Linux64RHEL6/TSCPBClassDict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableTSCPBClassDict();

extern "C" void G__set_cpp_environmentTSCPBClassDict() {
  G__cpp_reset_tagtableTSCPBClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTSCPBClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TSCPBClass */
static int G__TSCPBClassDict_183_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TSCPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TSCPBClass[n];
     } else {
       p = new((void*) gvp) TSCPBClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TSCPBClass;
     } else {
       p = new((void*) gvp) TSCPBClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TSCPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TSCPBClass((TSCPBClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TSCPBClass((TSCPBClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TSCPBClass*) G__getstructoffset())->GetScpdht());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TSCPBClass*) G__getstructoffset())->GetEdep());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TSCPBClass*) G__getstructoffset())->GetTime());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TSCPBClass*) G__getstructoffset())->GetPath());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TSCPBClass*) G__getstructoffset())->GetChi2());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TSCPBClass*) G__getstructoffset())->GetStatus());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TSCPBClass*) G__getstructoffset())->GetSector());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TSCPBClass*) G__getstructoffset())->GetPaddle());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TSCPBClass*) G__getstructoffset())->GetHit());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TSCPBClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TSCPBClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TSCPBClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TSCPBClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TSCPBClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TSCPBClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TSCPBClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TSCPBClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TSCPBClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TSCPBClassDict_183_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TSCPBClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TSCPBClassDict_183_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TSCPBClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TSCPBClass(*(TSCPBClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TSCPBClass G__TTSCPBClass;
static int G__TSCPBClassDict_183_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (TSCPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TSCPBClass*) (soff+(sizeof(TSCPBClass)*i)))->~G__TTSCPBClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TSCPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TSCPBClass*) (soff))->~G__TTSCPBClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TSCPBClassDict_183_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TSCPBClass* dest = (TSCPBClass*) G__getstructoffset();
   *dest = *(TSCPBClass*) libp->para[0].ref;
   const TSCPBClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TSCPBClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTSCPBClassDict {
 public:
  G__Sizep2memfuncTSCPBClassDict(): p(&G__Sizep2memfuncTSCPBClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTSCPBClassDict::*p)();
};

size_t G__get_sizep2memfuncTSCPBClassDict()
{
  G__Sizep2memfuncTSCPBClassDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceTSCPBClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass))) {
     TSCPBClass *G__Lderived;
     G__Lderived=(TSCPBClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass),G__get_linked_tagnum(&G__TSCPBClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTSCPBClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Float_t",102,-1,0,-1);
   G__setnewtype(-1,"Float 4 bytes (float)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TSCPBClass */
static void G__setup_memvarTSCPBClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass));
   { TSCPBClass *p; p=(TSCPBClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->Scpdht)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Scpdht=",0,"10000*sector+100*SC_PD_ID+Hit_ID in SCR ");
   G__memvar_setup((void*)((long)(&p->Edep)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Edep=",0,"Deposited energy (dE/dX)");
   G__memvar_setup((void*)((long)(&p->Time)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Time=",0,"Flight time relative to the evnt start time");
   G__memvar_setup((void*)((long)(&p->Path)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Path=",0,"Path lenght from target");
   G__memvar_setup((void*)((long)(&p->Chi2sc)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Chi2sc=",0,"Quality measure of geometrical matching");
   G__memvar_setup((void*)((long)(&p->Status)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Status=",0,"Status word (not defined yet)");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TSCPBClassDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTSCPBClassDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTSCPBClass(void) {
   /* TSCPBClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass));
   G__memfunc_setup("TSCPBClass",882,G__TSCPBClassDict_183_0_1, 105, G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TSCPBClass",882,G__TSCPBClassDict_183_0_2, 105, G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass), -1, 0, 1, 1, 1, 0, "U 'TSCPBClass' - 0 - TmpSCPB", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetScpdht",902,G__TSCPBClassDict_183_0_3, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetEdep",670,G__TSCPBClassDict_183_0_4, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetTime",687,G__TSCPBClassDict_183_0_5, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetPath",685,G__TSCPBClassDict_183_0_6, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetChi2",614,G__TSCPBClassDict_183_0_7, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetStatus",932,G__TSCPBClassDict_183_0_8, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetSector",912,G__TSCPBClassDict_183_0_9, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetPaddle",874,G__TSCPBClassDict_183_0_10, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHit",581,G__TSCPBClassDict_183_0_11, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__TSCPBClassDict_183_0_12, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TSCPBClassDict_183_0_13, 85, G__get_linked_tagnum(&G__TSCPBClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TSCPBClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TSCPBClassDict_183_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TSCPBClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TSCPBClassDict_183_0_15, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TSCPBClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TSCPBClassDict_183_0_16, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TSCPBClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TSCPBClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TSCPBClassDict_183_0_20, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TSCPBClassDict_183_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TSCPBClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TSCPBClassDict_183_0_22, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TSCPBClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TSCPBClassDict_183_0_23, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TSCPBClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TSCPBClassDict_183_0_24, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TSCPBClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TSCPBClass", 882, G__TSCPBClassDict_183_0_25, (int) ('i'), G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass), -1, 0, 1, 1, 1, 0, "u 'TSCPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TSCPBClass", 1008, G__TSCPBClassDict_183_0_26, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TSCPBClassDict_183_0_27, (int) ('u'), G__get_linked_tagnum(&G__TSCPBClassDictLN_TSCPBClass), -1, 1, 1, 1, 1, 0, "u 'TSCPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTSCPBClassDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalTSCPBClassDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTSCPBClassDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TSCPBClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TSCPBClassDictLN_TSCPBClass = { "TSCPBClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTSCPBClassDict() {
  G__TSCPBClassDictLN_TClass.tagnum = -1 ;
  G__TSCPBClassDictLN_TBuffer.tagnum = -1 ;
  G__TSCPBClassDictLN_TMemberInspector.tagnum = -1 ;
  G__TSCPBClassDictLN_TObject.tagnum = -1 ;
  G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TSCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TSCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TSCPBClassDictLN_TSCPBClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTSCPBClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TSCPBClassDictLN_TSCPBClass),sizeof(TSCPBClass),-1,62720,"Class for accessing the SCPB bank: Time of Flight.",G__setup_memvarTSCPBClass,G__setup_memfuncTSCPBClass);
}
extern "C" void G__cpp_setupTSCPBClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTSCPBClassDict()");
  G__set_cpp_environmentTSCPBClassDict();
  G__cpp_setup_tagtableTSCPBClassDict();

  G__cpp_setup_inheritanceTSCPBClassDict();

  G__cpp_setup_typetableTSCPBClassDict();

  G__cpp_setup_memvarTSCPBClassDict();

  G__cpp_setup_memfuncTSCPBClassDict();
  G__cpp_setup_globalTSCPBClassDict();
  G__cpp_setup_funcTSCPBClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTSCPBClassDict();
  return;
}
class G__cpp_setup_initTSCPBClassDict {
  public:
    G__cpp_setup_initTSCPBClassDict() { G__add_setup_func("TSCPBClassDict",(G__incsetup)(&G__cpp_setupTSCPBClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTSCPBClassDict() { G__remove_setup_func("TSCPBClassDict"); }
};
G__cpp_setup_initTSCPBClassDict G__cpp_setup_initializerTSCPBClassDict;

