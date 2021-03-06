//
// File generated by /u/apps/root/5.34.36/root/bin/rootcint at Fri Feb 24 17:38:29 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL7dITCCPBClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TCCPBClassDict.h"

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
   void TCCPBClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TCCPBClass(void *p = 0);
   static void *newArray_TCCPBClass(Long_t size, void *p);
   static void delete_TCCPBClass(void *p);
   static void deleteArray_TCCPBClass(void *p);
   static void destruct_TCCPBClass(void *p);
   static void streamer_TCCPBClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::TCCPBClass*)
   {
      ::TCCPBClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TCCPBClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TCCPBClass", ::TCCPBClass::Class_Version(), "./TCCPBClass.h", 22,
                  typeid(::TCCPBClass), ::ROOT::DefineBehavior(ptr, ptr),
                  &::TCCPBClass::Dictionary, isa_proxy, 0,
                  sizeof(::TCCPBClass) );
      instance.SetNew(&new_TCCPBClass);
      instance.SetNewArray(&newArray_TCCPBClass);
      instance.SetDelete(&delete_TCCPBClass);
      instance.SetDeleteArray(&deleteArray_TCCPBClass);
      instance.SetDestructor(&destruct_TCCPBClass);
      instance.SetStreamerFunc(&streamer_TCCPBClass);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::TCCPBClass*)
   {
      return GenerateInitInstanceLocal((::TCCPBClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCCPBClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr TCCPBClass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TCCPBClass::Class_Name()
{
   return "TCCPBClass";
}

//______________________________________________________________________________
const char *TCCPBClass::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TCCPBClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TCCPBClass::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TCCPBClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TCCPBClass::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TCCPBClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TCCPBClass::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TCCPBClass*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void TCCPBClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class TCCPBClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Scsght;
      R__b >> Nphe;
      R__b >> Time;
      R__b >> Path;
      R__b >> Chi2cc;
      R__b >> Status;
      R__b.CheckByteCount(R__s, R__c, TCCPBClass::IsA());
   } else {
      R__c = R__b.WriteVersion(TCCPBClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Scsght;
      R__b << Nphe;
      R__b << Time;
      R__b << Path;
      R__b << Chi2cc;
      R__b << Status;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TCCPBClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TCCPBClass.
      TClass *R__cl = ::TCCPBClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Scsght", &Scsght);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Nphe", &Nphe);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Time", &Time);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Path", &Path);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Chi2cc", &Chi2cc);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Status", &Status);
      TObject::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrappers around operator new
   static void *new_TCCPBClass(void *p) {
      return  p ? new(p) ::TCCPBClass : new ::TCCPBClass;
   }
   static void *newArray_TCCPBClass(Long_t nElements, void *p) {
      return p ? new(p) ::TCCPBClass[nElements] : new ::TCCPBClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCCPBClass(void *p) {
      delete ((::TCCPBClass*)p);
   }
   static void deleteArray_TCCPBClass(void *p) {
      delete [] ((::TCCPBClass*)p);
   }
   static void destruct_TCCPBClass(void *p) {
      typedef ::TCCPBClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TCCPBClass(TBuffer &buf, void *obj) {
      ((::TCCPBClass*)obj)->::TCCPBClass::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::TCCPBClass

/********************************************************
* dicts/Linux64RHEL7/TCCPBClassDict.cc
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

extern "C" void G__cpp_reset_tagtableTCCPBClassDict();

extern "C" void G__set_cpp_environmentTCCPBClassDict() {
  G__cpp_reset_tagtableTCCPBClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTCCPBClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TCCPBClass */
static int G__TCCPBClassDict_184_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TCCPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TCCPBClass[n];
     } else {
       p = new((void*) gvp) TCCPBClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TCCPBClass;
     } else {
       p = new((void*) gvp) TCCPBClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TCCPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TCCPBClass((TCCPBClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TCCPBClass((TCCPBClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TCCPBClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TCCPBClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TCCPBClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TCCPBClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TCCPBClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TCCPBClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TCCPBClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TCCPBClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TCCPBClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TCCPBClassDict_184_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TCCPBClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TCCPBClassDict_184_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TCCPBClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TCCPBClass(*(TCCPBClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TCCPBClass G__TTCCPBClass;
static int G__TCCPBClassDict_184_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (TCCPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TCCPBClass*) (soff+(sizeof(TCCPBClass)*i)))->~G__TTCCPBClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TCCPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TCCPBClass*) (soff))->~G__TTCCPBClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TCCPBClassDict_184_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TCCPBClass* dest = (TCCPBClass*) G__getstructoffset();
   *dest = *(TCCPBClass*) libp->para[0].ref;
   const TCCPBClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TCCPBClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTCCPBClassDict {
 public:
  G__Sizep2memfuncTCCPBClassDict(): p(&G__Sizep2memfuncTCCPBClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTCCPBClassDict::*p)();
};

size_t G__get_sizep2memfuncTCCPBClassDict()
{
  G__Sizep2memfuncTCCPBClassDict a;
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
extern "C" void G__cpp_setup_inheritanceTCCPBClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass))) {
     TCCPBClass *G__Lderived;
     G__Lderived=(TCCPBClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass),G__get_linked_tagnum(&G__TCCPBClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTCCPBClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TCCPBClass */
static void G__setup_memvarTCCPBClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass));
   { TCCPBClass *p; p=(TCCPBClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->Scsght)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Scsght=",0,"1000*sector+100*CC_segm_ID+Hit_ID in CCRC");
   G__memvar_setup((void*)((long)(&p->Nphe)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Nphe=",0,"Number of photo-electrons");
   G__memvar_setup((void*)((long)(&p->Time)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Time=",0,"Flight time relative to the evnt start time");
   G__memvar_setup((void*)((long)(&p->Path)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Path=",0,"Path lenght from target");
   G__memvar_setup((void*)((long)(&p->Chi2cc)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Chi2cc=",0,"Quality measure of geometrical matching");
   G__memvar_setup((void*)((long)(&p->Status)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"Status=",0,"Status word (not defined yet)");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TCCPBClassDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTCCPBClassDict() {
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
static void G__setup_memfuncTCCPBClass(void) {
   /* TCCPBClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass));
   G__memfunc_setup("TCCPBClass",866,G__TCCPBClassDict_184_0_1, 105, G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TCCPBClass",866,G__TCCPBClassDict_184_0_2, 105, G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass), -1, 0, 1, 1, 1, 0, "U 'TCCPBClass' - 0 - TmpCCPB", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__TCCPBClassDict_184_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TCCPBClassDict_184_0_4, 85, G__get_linked_tagnum(&G__TCCPBClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TCCPBClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TCCPBClassDict_184_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TCCPBClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TCCPBClassDict_184_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TCCPBClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TCCPBClassDict_184_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TCCPBClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TCCPBClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TCCPBClassDict_184_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TCCPBClassDict_184_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TCCPBClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TCCPBClassDict_184_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TCCPBClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TCCPBClassDict_184_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TCCPBClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TCCPBClassDict_184_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TCCPBClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TCCPBClass", 866, G__TCCPBClassDict_184_0_16, (int) ('i'), G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass), -1, 0, 1, 1, 1, 0, "u 'TCCPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TCCPBClass", 992, G__TCCPBClassDict_184_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TCCPBClassDict_184_0_18, (int) ('u'), G__get_linked_tagnum(&G__TCCPBClassDictLN_TCCPBClass), -1, 1, 1, 1, 1, 0, "u 'TCCPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTCCPBClassDict() {
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
extern "C" void G__cpp_setup_globalTCCPBClassDict() {
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

extern "C" void G__cpp_setup_funcTCCPBClassDict() {
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
G__linked_taginfo G__TCCPBClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TCCPBClassDictLN_TCCPBClass = { "TCCPBClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTCCPBClassDict() {
  G__TCCPBClassDictLN_TClass.tagnum = -1 ;
  G__TCCPBClassDictLN_TBuffer.tagnum = -1 ;
  G__TCCPBClassDictLN_TMemberInspector.tagnum = -1 ;
  G__TCCPBClassDictLN_TObject.tagnum = -1 ;
  G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TCCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TCCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TCCPBClassDictLN_TCCPBClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTCCPBClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TCCPBClassDictLN_TCCPBClass),sizeof(TCCPBClass),-1,62720,"Class for accessing the CCPB bank: Cherencov.",G__setup_memvarTCCPBClass,G__setup_memfuncTCCPBClass);
}
extern "C" void G__cpp_setupTCCPBClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTCCPBClassDict()");
  G__set_cpp_environmentTCCPBClassDict();
  G__cpp_setup_tagtableTCCPBClassDict();

  G__cpp_setup_inheritanceTCCPBClassDict();

  G__cpp_setup_typetableTCCPBClassDict();

  G__cpp_setup_memvarTCCPBClassDict();

  G__cpp_setup_memfuncTCCPBClassDict();
  G__cpp_setup_globalTCCPBClassDict();
  G__cpp_setup_funcTCCPBClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTCCPBClassDict();
  return;
}
class G__cpp_setup_initTCCPBClassDict {
  public:
    G__cpp_setup_initTCCPBClassDict() { G__add_setup_func("TCCPBClassDict",(G__incsetup)(&G__cpp_setupTCCPBClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTCCPBClassDict() { G__remove_setup_func("TCCPBClassDict"); }
};
G__cpp_setup_initTCCPBClassDict G__cpp_setup_initializerTCCPBClassDict;

