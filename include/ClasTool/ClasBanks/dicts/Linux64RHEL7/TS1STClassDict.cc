//
// File generated by /u/apps/root/5.34.36/root/bin/rootcint at Fri Feb 24 17:38:16 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL7dITS1STClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TS1STClassDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
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

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

namespace ROOTDict {
   void TS1STClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TS1STClass(void *p = 0);
   static void *newArray_TS1STClass(Long_t size, void *p);
   static void delete_TS1STClass(void *p);
   static void deleteArray_TS1STClass(void *p);
   static void destruct_TS1STClass(void *p);
   static void streamer_TS1STClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::TS1STClass*)
   {
      ::TS1STClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TS1STClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TS1STClass", ::TS1STClass::Class_Version(), "./TS1STClass.h", 20,
                  typeid(::TS1STClass), ::ROOT::DefineBehavior(ptr, ptr),
                  &::TS1STClass::Dictionary, isa_proxy, 0,
                  sizeof(::TS1STClass) );
      instance.SetNew(&new_TS1STClass);
      instance.SetNewArray(&newArray_TS1STClass);
      instance.SetDelete(&delete_TS1STClass);
      instance.SetDeleteArray(&deleteArray_TS1STClass);
      instance.SetDestructor(&destruct_TS1STClass);
      instance.SetStreamerFunc(&streamer_TS1STClass);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::TS1STClass*)
   {
      return GenerateInitInstanceLocal((::TS1STClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TS1STClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr TS1STClass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TS1STClass::Class_Name()
{
   return "TS1STClass";
}

//______________________________________________________________________________
const char *TS1STClass::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TS1STClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TS1STClass::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TS1STClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TS1STClass::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TS1STClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TS1STClass::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TS1STClass*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void TS1STClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class TS1STClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Latch1_bit1_count;
      R__b >> Latch1_bit2_count;
      R__b >> Latch1_bit3_count;
      R__b >> Latch1_bit4_count;
      R__b >> Latch1_bit5_count;
      R__b >> Latch1_bit6_count;
      R__b >> Latch1_bit7_count;
      R__b >> Latch1_bit8_count;
      R__b >> Latch1_bit9_count;
      R__b >> Latch1_bit10_count;
      R__b >> Latch1_bit11_count;
      R__b >> Latch1_bit12_count;
      R__b >> Event_count;
      R__b >> Latch1_zero_count;
      R__b >> Latch1_empty_count;
      R__b >> Latch1_not_empty_count;
      R__b >> Latch1_ok_count;
      R__b >> Level2_sector1;
      R__b >> Level2_sector2;
      R__b >> Level2_sector3;
      R__b >> Level2_sector4;
      R__b >> Level2_sector5;
      R__b >> Level2_sector6;
      R__b >> Level2_pass;
      R__b >> Level2_fail;
      R__b >> Latch2_zero_count;
      R__b >> Latch2_empty_count;
      R__b >> Latch2_not_empty_count;
      R__b >> Latch2_ok_count;
      R__b >> Roc_13_count;
      R__b >> Roc_15_count;
      R__b >> L1l2_zero_count;
      R__b >> L1zero_13_count;
      R__b >> L2zero_13_count;
      R__b >> L1l2zero_13_count;
      R__b.CheckByteCount(R__s, R__c, TS1STClass::IsA());
   } else {
      R__c = R__b.WriteVersion(TS1STClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Latch1_bit1_count;
      R__b << Latch1_bit2_count;
      R__b << Latch1_bit3_count;
      R__b << Latch1_bit4_count;
      R__b << Latch1_bit5_count;
      R__b << Latch1_bit6_count;
      R__b << Latch1_bit7_count;
      R__b << Latch1_bit8_count;
      R__b << Latch1_bit9_count;
      R__b << Latch1_bit10_count;
      R__b << Latch1_bit11_count;
      R__b << Latch1_bit12_count;
      R__b << Event_count;
      R__b << Latch1_zero_count;
      R__b << Latch1_empty_count;
      R__b << Latch1_not_empty_count;
      R__b << Latch1_ok_count;
      R__b << Level2_sector1;
      R__b << Level2_sector2;
      R__b << Level2_sector3;
      R__b << Level2_sector4;
      R__b << Level2_sector5;
      R__b << Level2_sector6;
      R__b << Level2_pass;
      R__b << Level2_fail;
      R__b << Latch2_zero_count;
      R__b << Latch2_empty_count;
      R__b << Latch2_not_empty_count;
      R__b << Latch2_ok_count;
      R__b << Roc_13_count;
      R__b << Roc_15_count;
      R__b << L1l2_zero_count;
      R__b << L1zero_13_count;
      R__b << L2zero_13_count;
      R__b << L1l2zero_13_count;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TS1STClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TS1STClass.
      TClass *R__cl = ::TS1STClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit1_count", &Latch1_bit1_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit2_count", &Latch1_bit2_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit3_count", &Latch1_bit3_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit4_count", &Latch1_bit4_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit5_count", &Latch1_bit5_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit6_count", &Latch1_bit6_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit7_count", &Latch1_bit7_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit8_count", &Latch1_bit8_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit9_count", &Latch1_bit9_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit10_count", &Latch1_bit10_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit11_count", &Latch1_bit11_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_bit12_count", &Latch1_bit12_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Event_count", &Event_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_zero_count", &Latch1_zero_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_empty_count", &Latch1_empty_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_not_empty_count", &Latch1_not_empty_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch1_ok_count", &Latch1_ok_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector1", &Level2_sector1);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector2", &Level2_sector2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector3", &Level2_sector3);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector4", &Level2_sector4);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector5", &Level2_sector5);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_sector6", &Level2_sector6);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_pass", &Level2_pass);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Level2_fail", &Level2_fail);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch2_zero_count", &Latch2_zero_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch2_empty_count", &Latch2_empty_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch2_not_empty_count", &Latch2_not_empty_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Latch2_ok_count", &Latch2_ok_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Roc_13_count", &Roc_13_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Roc_15_count", &Roc_15_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "L1l2_zero_count", &L1l2_zero_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "L1zero_13_count", &L1zero_13_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "L2zero_13_count", &L2zero_13_count);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "L1l2zero_13_count", &L1l2zero_13_count);
      TObject::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrappers around operator new
   static void *new_TS1STClass(void *p) {
      return  p ? new(p) ::TS1STClass : new ::TS1STClass;
   }
   static void *newArray_TS1STClass(Long_t nElements, void *p) {
      return p ? new(p) ::TS1STClass[nElements] : new ::TS1STClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_TS1STClass(void *p) {
      delete ((::TS1STClass*)p);
   }
   static void deleteArray_TS1STClass(void *p) {
      delete [] ((::TS1STClass*)p);
   }
   static void destruct_TS1STClass(void *p) {
      typedef ::TS1STClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TS1STClass(TBuffer &buf, void *obj) {
      ((::TS1STClass*)obj)->::TS1STClass::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::TS1STClass

/********************************************************
* dicts/Linux64RHEL7/TS1STClassDict.cc
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

extern "C" void G__cpp_reset_tagtableTS1STClassDict();

extern "C" void G__set_cpp_environmentTS1STClassDict() {
  G__cpp_reset_tagtableTS1STClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTS1STClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TS1STClass */
static int G__TS1STClassDict_184_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TS1STClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TS1STClass[n];
     } else {
       p = new((void*) gvp) TS1STClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TS1STClass;
     } else {
       p = new((void*) gvp) TS1STClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TS1STClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TS1STClass((TS1STClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TS1STClass((TS1STClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TS1STClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TS1STClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TS1STClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TS1STClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TS1STClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TS1STClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TS1STClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TS1STClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TS1STClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TS1STClassDict_184_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TS1STClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TS1STClassDict_184_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TS1STClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TS1STClass(*(TS1STClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TS1STClass G__TTS1STClass;
static int G__TS1STClassDict_184_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (TS1STClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TS1STClass*) (soff+(sizeof(TS1STClass)*i)))->~G__TTS1STClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TS1STClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TS1STClass*) (soff))->~G__TTS1STClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TS1STClassDict_184_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TS1STClass* dest = (TS1STClass*) G__getstructoffset();
   *dest = *(TS1STClass*) libp->para[0].ref;
   const TS1STClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TS1STClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTS1STClassDict {
 public:
  G__Sizep2memfuncTS1STClassDict(): p(&G__Sizep2memfuncTS1STClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTS1STClassDict::*p)();
};

size_t G__get_sizep2memfuncTS1STClassDict()
{
  G__Sizep2memfuncTS1STClassDict a;
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
extern "C" void G__cpp_setup_inheritanceTS1STClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass))) {
     TS1STClass *G__Lderived;
     G__Lderived=(TS1STClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass),G__get_linked_tagnum(&G__TS1STClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTS1STClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TS1STClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TS1STClass */
static void G__setup_memvarTS1STClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass));
   { TS1STClass *p; p=(TS1STClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->Latch1_bit1_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit1_count=",0,"Count trigger bit 1  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit2_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit2_count=",0,"Count trigger bit 2  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit3_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit3_count=",0,"Count trigger bit 3  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit4_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit4_count=",0,"Count trigger bit 4  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit5_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit5_count=",0,"Count trigger bit 5  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit6_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit6_count=",0,"Count trigger bit 6  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit7_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit7_count=",0,"Count trigger bit 7  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit8_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit8_count=",0,"Count trigger bit 8  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit9_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit9_count=",0,"Count trigger bit 9  latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit10_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit10_count=",0,"Count trigger bit 10 latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit11_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit11_count=",0,"Count trigger bit 11 latched events            ");
   G__memvar_setup((void*)((long)(&p->Latch1_bit12_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_bit12_count=",0,"Count trigger bit 12 latched events            ");
   G__memvar_setup((void*)((long)(&p->Event_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Event_count=",0,"Latched event count this run                   ");
   G__memvar_setup((void*)((long)(&p->Latch1_zero_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_zero_count=",0,"Latch1 zero count (illegal)                    ");
   G__memvar_setup((void*)((long)(&p->Latch1_empty_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_empty_count=",0,"Latch1 empty count (illegal)                   ");
   G__memvar_setup((void*)((long)(&p->Latch1_not_empty_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_not_empty_count=",0,"Latch1 not empty on sync event (illegal)       ");
   G__memvar_setup((void*)((long)(&p->Latch1_ok_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch1_ok_count=",0,"Latch1 ok                                      ");
   G__memvar_setup((void*)((long)(&p->Level2_sector1)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector1=",0,"Level2 sector1 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_sector2)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector2=",0,"Level2 sector2 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_sector3)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector3=",0,"Level2 sector3 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_sector4)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector4=",0,"Level2 sector4 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_sector5)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector5=",0,"Level2 sector5 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_sector6)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_sector6=",0,"Level2 sector6 count                           ");
   G__memvar_setup((void*)((long)(&p->Level2_pass)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_pass=",0,"Level2 pass count                              ");
   G__memvar_setup((void*)((long)(&p->Level2_fail)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Level2_fail=",0,"Level2 fail count                              ");
   G__memvar_setup((void*)((long)(&p->Latch2_zero_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch2_zero_count=",0,"Latch2 zero count (illegal)                    ");
   G__memvar_setup((void*)((long)(&p->Latch2_empty_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch2_empty_count=",0,"Latch2 empty count (illegal)                   ");
   G__memvar_setup((void*)((long)(&p->Latch2_not_empty_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch2_not_empty_count=",0,"Latch2 not empty on sync event (illegal)       ");
   G__memvar_setup((void*)((long)(&p->Latch2_ok_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Latch2_ok_count=",0,"Latch2 ok                                      ");
   G__memvar_setup((void*)((long)(&p->Roc_13_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Roc_13_count=",0,"Roc code 13 count (zero latch)                 ");
   G__memvar_setup((void*)((long)(&p->Roc_15_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Roc_15_count=",0,"Roc code 15 count (illegal)                    ");
   G__memvar_setup((void*)((long)(&p->L1l2_zero_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"L1l2_zero_count=",0,"(latch1==0)&&(latch2==0)                       ");
   G__memvar_setup((void*)((long)(&p->L1zero_13_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"L1zero_13_count=",0,"(latch1==0)&&(roc_code==13)                    ");
   G__memvar_setup((void*)((long)(&p->L2zero_13_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"L2zero_13_count=",0,"(latch2==0)&&(roc_code==13)                    ");
   G__memvar_setup((void*)((long)(&p->L1l2zero_13_count)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"L1l2zero_13_count=",0,"(latch1==0)&&(latch2==0)&&(roc_code==13)       ");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TS1STClassDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTS1STClassDict() {
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
static void G__setup_memfuncTS1STClass(void) {
   /* TS1STClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass));
   G__memfunc_setup("TS1STClass",885,G__TS1STClassDict_184_0_1, 105, G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TS1STClass",885,G__TS1STClassDict_184_0_2, 105, G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass), -1, 0, 1, 1, 1, 0, "U 'TS1STClass' - 0 - TmpS1ST", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__TS1STClassDict_184_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TS1STClassDict_184_0_4, 85, G__get_linked_tagnum(&G__TS1STClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TS1STClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TS1STClassDict_184_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TS1STClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TS1STClassDict_184_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TS1STClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TS1STClassDict_184_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TS1STClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TS1STClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TS1STClassDict_184_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TS1STClassDict_184_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TS1STClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TS1STClassDict_184_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TS1STClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TS1STClassDict_184_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TS1STClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TS1STClassDict_184_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TS1STClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TS1STClass", 885, G__TS1STClassDict_184_0_16, (int) ('i'), G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass), -1, 0, 1, 1, 1, 0, "u 'TS1STClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TS1STClass", 1011, G__TS1STClassDict_184_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TS1STClassDict_184_0_18, (int) ('u'), G__get_linked_tagnum(&G__TS1STClassDictLN_TS1STClass), -1, 1, 1, 1, 1, 0, "u 'TS1STClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTS1STClassDict() {
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
extern "C" void G__cpp_setup_globalTS1STClassDict() {
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
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTS1STClassDict() {
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
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TS1STClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TS1STClassDictLN_TS1STClass = { "TS1STClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTS1STClassDict() {
  G__TS1STClassDictLN_TClass.tagnum = -1 ;
  G__TS1STClassDictLN_TBuffer.tagnum = -1 ;
  G__TS1STClassDictLN_TMemberInspector.tagnum = -1 ;
  G__TS1STClassDictLN_TObject.tagnum = -1 ;
  G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TS1STClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TS1STClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TS1STClassDictLN_TS1STClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTS1STClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TS1STClassDictLN_TS1STClass),sizeof(TS1STClass),-1,62720,"Scaler bank with trigger and latch information",G__setup_memvarTS1STClass,G__setup_memfuncTS1STClass);
}
extern "C" void G__cpp_setupTS1STClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTS1STClassDict()");
  G__set_cpp_environmentTS1STClassDict();
  G__cpp_setup_tagtableTS1STClassDict();

  G__cpp_setup_inheritanceTS1STClassDict();

  G__cpp_setup_typetableTS1STClassDict();

  G__cpp_setup_memvarTS1STClassDict();

  G__cpp_setup_memfuncTS1STClassDict();
  G__cpp_setup_globalTS1STClassDict();
  G__cpp_setup_funcTS1STClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTS1STClassDict();
  return;
}
class G__cpp_setup_initTS1STClassDict {
  public:
    G__cpp_setup_initTS1STClassDict() { G__add_setup_func("TS1STClassDict",(G__incsetup)(&G__cpp_setupTS1STClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTS1STClassDict() { G__remove_setup_func("TS1STClassDict"); }
};
G__cpp_setup_initTS1STClassDict G__cpp_setup_initializerTS1STClassDict;

