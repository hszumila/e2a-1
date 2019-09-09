//
// File generated by /u/apps/root/5.34.13/root/bin/rootcint at Fri Jan 20 19:15:29 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL6dITHEADClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "THEADClassDict.h"

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
   void THEADClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_THEADClass(void *p = 0);
   static void *newArray_THEADClass(Long_t size, void *p);
   static void delete_THEADClass(void *p);
   static void deleteArray_THEADClass(void *p);
   static void destruct_THEADClass(void *p);
   static void streamer_THEADClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::THEADClass*)
   {
      ::THEADClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::THEADClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("THEADClass", ::THEADClass::Class_Version(), "./THEADClass.h", 17,
                  typeid(::THEADClass), DefineBehavior(ptr, ptr),
                  &::THEADClass::Dictionary, isa_proxy, 0,
                  sizeof(::THEADClass) );
      instance.SetNew(&new_THEADClass);
      instance.SetNewArray(&newArray_THEADClass);
      instance.SetDelete(&delete_THEADClass);
      instance.SetDeleteArray(&deleteArray_THEADClass);
      instance.SetDestructor(&destruct_THEADClass);
      instance.SetStreamerFunc(&streamer_THEADClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::THEADClass*)
   {
      return GenerateInitInstanceLocal((::THEADClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::THEADClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *THEADClass::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *THEADClass::Class_Name()
{
   return "THEADClass";
}

//______________________________________________________________________________
const char *THEADClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THEADClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int THEADClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THEADClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void THEADClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THEADClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *THEADClass::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THEADClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void THEADClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class THEADClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Version;
      R__b >> Nrun;
      R__b >> Nevent;
      R__b >> Time;
      R__b >> Type;
      R__b >> Roc;
      R__b >> Evtclass;
      R__b >> Trigbits;
      R__b.CheckByteCount(R__s, R__c, THEADClass::IsA());
   } else {
      R__c = R__b.WriteVersion(THEADClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Version;
      R__b << Nrun;
      R__b << Nevent;
      R__b << Time;
      R__b << Type;
      R__b << Roc;
      R__b << Evtclass;
      R__b << Trigbits;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void THEADClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class THEADClass.
      TClass *R__cl = ::THEADClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Version", &Version);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Nrun", &Nrun);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Nevent", &Nevent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Time", &Time);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Type", &Type);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Roc", &Roc);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Evtclass", &Evtclass);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Trigbits", &Trigbits);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_THEADClass(void *p) {
      return  p ? new(p) ::THEADClass : new ::THEADClass;
   }
   static void *newArray_THEADClass(Long_t nElements, void *p) {
      return p ? new(p) ::THEADClass[nElements] : new ::THEADClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_THEADClass(void *p) {
      delete ((::THEADClass*)p);
   }
   static void deleteArray_THEADClass(void *p) {
      delete [] ((::THEADClass*)p);
   }
   static void destruct_THEADClass(void *p) {
      typedef ::THEADClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_THEADClass(TBuffer &buf, void *obj) {
      ((::THEADClass*)obj)->::THEADClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::THEADClass

/********************************************************
* dicts/Linux64RHEL6/THEADClassDict.cc
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

extern "C" void G__cpp_reset_tagtableTHEADClassDict();

extern "C" void G__set_cpp_environmentTHEADClassDict() {
  G__cpp_reset_tagtableTHEADClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTHEADClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* THEADClass */
static int G__THEADClassDict_183_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THEADClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new THEADClass[n];
     } else {
       p = new((void*) gvp) THEADClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new THEADClass;
     } else {
       p = new((void*) gvp) THEADClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THEADClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new THEADClass((THEADClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) THEADClass((THEADClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((THEADClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) THEADClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THEADClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) THEADClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      THEADClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((THEADClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THEADClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) THEADClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THEADClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THEADClassDict_183_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) THEADClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__THEADClassDict_183_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   THEADClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new THEADClass(*(THEADClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef THEADClass G__TTHEADClass;
static int G__THEADClassDict_183_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (THEADClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((THEADClass*) (soff+(sizeof(THEADClass)*i)))->~G__TTHEADClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (THEADClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((THEADClass*) (soff))->~G__TTHEADClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__THEADClassDict_183_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THEADClass* dest = (THEADClass*) G__getstructoffset();
   *dest = *(THEADClass*) libp->para[0].ref;
   const THEADClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* THEADClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTHEADClassDict {
 public:
  G__Sizep2memfuncTHEADClassDict(): p(&G__Sizep2memfuncTHEADClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTHEADClassDict::*p)();
};

size_t G__get_sizep2memfuncTHEADClassDict()
{
  G__Sizep2memfuncTHEADClassDict a;
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
extern "C" void G__cpp_setup_inheritanceTHEADClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass))) {
     THEADClass *G__Lderived;
     G__Lderived=(THEADClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass),G__get_linked_tagnum(&G__THEADClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTHEADClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__THEADClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__THEADClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__THEADClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__THEADClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* THEADClass */
static void G__setup_memvarTHEADClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass));
   { THEADClass *p; p=(THEADClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->Version)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Version=",0,"Version number from Recsis");
   G__memvar_setup((void*)((long)(&p->Nrun)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Nrun=",0,"Run number for this run");
   G__memvar_setup((void*)((long)(&p->Nevent)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Nevent=",0,"Event number for this run, starts with 1, increases with one for BOS each event in BOS file.");
   G__memvar_setup((void*)((long)(&p->Time)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Time=",0,"Event Time as UNIX time.");
   G__memvar_setup((void*)((long)(&p->Type)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Type=",0,"Event Type: 1-9 Physics event (2= sync, 4=level2 late fail) 10 Scaler event. < 0 Monte Carlo");
   G__memvar_setup((void*)((long)(&p->Roc)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Roc=",0,"=0 Sync status ok, >0 = bit pattern of offending ROC.");
   G__memvar_setup((void*)((long)(&p->Evtclass)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Evtclass=",0,"Event Classification from DAQ 0=special event, 1-15 Physics event, 16 Sync Event, 17 Prestart, 18 Go, 19 Pause, 20 End.");
   G__memvar_setup((void*)((long)(&p->Trigbits)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Trigbits=",0,"Level 1 Trigger Latch word.");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__THEADClassDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTHEADClassDict() {
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
static void G__setup_memfuncTHEADClass(void) {
   /* THEADClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass));
   G__memfunc_setup("THEADClass",860,G__THEADClassDict_183_0_1, 105, G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("THEADClass",860,G__THEADClassDict_183_0_2, 105, G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass), -1, 0, 1, 1, 1, 0, "U 'THEADClass' - 0 - TmpHEAD", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__THEADClassDict_183_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__THEADClassDict_183_0_4, 85, G__get_linked_tagnum(&G__THEADClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&THEADClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__THEADClassDict_183_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THEADClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__THEADClassDict_183_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&THEADClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__THEADClassDict_183_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&THEADClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__THEADClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__THEADClassDict_183_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__THEADClassDict_183_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THEADClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__THEADClassDict_183_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&THEADClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__THEADClassDict_183_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THEADClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__THEADClassDict_183_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&THEADClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("THEADClass", 860, G__THEADClassDict_183_0_16, (int) ('i'), G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass), -1, 0, 1, 1, 1, 0, "u 'THEADClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~THEADClass", 986, G__THEADClassDict_183_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__THEADClassDict_183_0_18, (int) ('u'), G__get_linked_tagnum(&G__THEADClassDictLN_THEADClass), -1, 1, 1, 1, 1, 0, "u 'THEADClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTHEADClassDict() {
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
extern "C" void G__cpp_setup_globalTHEADClassDict() {
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

extern "C" void G__cpp_setup_funcTHEADClassDict() {
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
G__linked_taginfo G__THEADClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__THEADClassDictLN_THEADClass = { "THEADClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTHEADClassDict() {
  G__THEADClassDictLN_TClass.tagnum = -1 ;
  G__THEADClassDictLN_TBuffer.tagnum = -1 ;
  G__THEADClassDictLN_TMemberInspector.tagnum = -1 ;
  G__THEADClassDictLN_TObject.tagnum = -1 ;
  G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__THEADClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__THEADClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__THEADClassDictLN_THEADClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTHEADClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__THEADClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__THEADClassDictLN_THEADClass),sizeof(THEADClass),-1,62720,"Header information for event",G__setup_memvarTHEADClass,G__setup_memfuncTHEADClass);
}
extern "C" void G__cpp_setupTHEADClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTHEADClassDict()");
  G__set_cpp_environmentTHEADClassDict();
  G__cpp_setup_tagtableTHEADClassDict();

  G__cpp_setup_inheritanceTHEADClassDict();

  G__cpp_setup_typetableTHEADClassDict();

  G__cpp_setup_memvarTHEADClassDict();

  G__cpp_setup_memfuncTHEADClassDict();
  G__cpp_setup_globalTHEADClassDict();
  G__cpp_setup_funcTHEADClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTHEADClassDict();
  return;
}
class G__cpp_setup_initTHEADClassDict {
  public:
    G__cpp_setup_initTHEADClassDict() { G__add_setup_func("THEADClassDict",(G__incsetup)(&G__cpp_setupTHEADClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTHEADClassDict() { G__remove_setup_func("THEADClassDict"); }
};
G__cpp_setup_initTHEADClassDict G__cpp_setup_initializerTHEADClassDict;
