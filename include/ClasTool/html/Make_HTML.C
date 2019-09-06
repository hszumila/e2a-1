{
// Make_HTML.C
//
// Author: Maurik Holtrop (UNH) 4/28/2000
//
// This is the master document generation script.
//
// This script contains a table of the directories where it will look for a
// file called HTML_Load.C. This file must contain instructions on how to 
// load the library of that directory, and which classes it contains, so that
// an html file can be made for each class.
//
// Look in E2AnaTool for an example HTML_Load.C script.
//

gROOT->Reset();

Int_t retval;
TString SourceDirs(".");

if(!gSystem->AccessPathName("../macros/Helper_Load_Libs.h",kReadPermission)){
  gROOT->LoadMacro("../macros/Helper_Load_Libs.h");
}else{
  cerr << "Could not load helper functions from ../macros/Helper_Load_Libs.h -- ABORT \n" ;
  gSystem->Exit(1);
}

if(gSystem->DynamicPathName("libPhysics.so",0) == 0){
  cerr << "Could not find the libPhysics.so library, which should be part or ROOT.\n" ;
  cerr << "Please check your ROOT installation, and make sure LD_LIBRARY_PATH is\n";
  cerr << "set correctly.\n\n";
  exit(1);
}else{
  if(!(retval=gSystem->Load(gSystem->DynamicPathName("libPhysics.so",0))))
    cerr << "Library libPhysics.so loaded." << endl;
  else{
    cerr << "Error loading libPhysics.so or already loaded: " << retval << endl;
  }
}

if(gSystem->DynamicPathName("libHtml.so",0) == 0){
  cerr << "Could not find the libHtml.so library, which should be part or ROOT.\n" ;
  cerr << "Please check your ROOT installation, and make sure LD_LIBRARY_PATH is\n";
  cerr << "set correctly.\n\n";
  exit(1);
}else{
  if(!(retval=gSystem->Load(gSystem->DynamicPathName("libHtml.so",0))))
    cerr << "Library libHtml.so loaded." << endl;
  else{
    cerr << "Error loading libHtml.so or already loaded: " << retval << endl;
  }
}

//
// This is where we generate the html documents.
//
THtml *mHtml= new THtml(); // Initialize the THtml class.

 mHtml->SetOutputDir("../html");
 mHtml->SetSourceDir(".:../ClasBanks:../VirtualReader:../DSTReader:../MapUtils:../TClasTool:../Utils:../PartSieve::../Maputils:../Filter:../NT10Reader:../Gsim");
////////////////////////////////////////////////////////////////////////////////
// 
// SECTION FOR INCLUDING DIRECTORIES INTO DOCUMENTATION BUILD.
//
// ADD NEW DIRECTORIES AT THE BOTTOM OF THE FOLLOWING LIST.
//
////////////////////////////////////////////////////////////////////////////////
//
// Load the class and documents, as instructed in the HTML_Load.C for 
// each directory. Make a documentation for each class.
//
// NOTE: If a class depends on another, it must come AFTER, so that the appropriate
// libs are loaded.
//
if(Load_A_Lib("libClasBanks.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libVirtualReader.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libDSTReader.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libClasTool.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libPartSieve.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libMapUtils.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(Load_A_Lib("libFilter.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
if(gROOT->GetClass("THbookFile")){
   if(Load_A_Lib("libNT10Reader.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
 }
TCCPBClass *ccpb = new TCCPBClass();       mHtml->MakeClass("TCCPBClass"); 
TDCPBClass *dcpb = new TDCPBClass();       mHtml->MakeClass("TDCPBClass"); 
TECPBClass *ecpb = new TECPBClass();       mHtml->MakeClass("TECPBClass"); 
TVirtualData *vdat = new TVirtualData();   mHtml->MakeClass("TVirtualData"); 
TEVNTClass *evnt = new TEVNTClass();       mHtml->MakeClass("TEVNTClass"); 
TGSIMClass *gsim = new TGSIMClass();       mHtml->MakeClass("TGSIMClass"); 
THEADClass *head = new THEADClass();       mHtml->MakeClass("THEADClass"); 
THEADERClass *header = new THEADERClass(); mHtml->MakeClass("THEADERClass"); 
THLSClass    *hls    = new THLSClass();    mHtml->MakeClass("THLSClass"); 
TLCPBClass *lcpb = new TLCPBClass();       mHtml->MakeClass("TLCPBClass"); 
TS1STClass *s1st = new TS1STClass();       mHtml->MakeClass("TS1STClass"); 
TSCPBClass *scpb = new TSCPBClass();       mHtml->MakeClass("TSCPBClass"); 
TTRGSClass *trgs = new TTRGSClass();       mHtml->MakeClass("TTRGSClass"); 
TTGBIClass *tgbi = new TTGBIClass();       mHtml->MakeClass("TTGBIClass"); 
TSTPBClass *stpb = new TSTPBClass();       mHtml->MakeClass("TSTPBClass"); 
TTGPBClass *tgpb = new TTGPBClass();       mHtml->MakeClass("TTGPBClass"); 
TTAGRClass *tagr = new TTAGRClass();       mHtml->MakeClass("TTAGRClass"); 
TVERTClass *vert = new TVERTClass();       mHtml->MakeClass("TVERTClass"); 
TMVRTClass *mvrt = new TMVRTClass();       mHtml->MakeClass("TMVRTClass"); 
TEPICClass *epic = new TEPICClass();       mHtml->MakeClass("TEPICClass"); 
TEPICValue *epval = new TEPICValue();      mHtml->MakeClass("TEPICValue"); 
TTBERClass *tber = new TTBERClass();       mHtml->MakeClass("TTBERClass");
TMCHDClass *mchd = new TMCHDClass();       mHtml->MakeClass("TMCHDClass"); 

TVirtualReader *virtread= new TVirtualReader(); mHtml->MakeClass("TVirtualReader");
 
TDSTReader *dstr = new TDSTReader();        mHtml->MakeClass("TDSTReader"); 
TMiniDSTReader *mdstr= new TMiniDSTReader(); mHtml->MakeClass("TMiniDSTReader"); 
TMiniDSTWriter *mdstw= new TMiniDSTWriter(); mHtml->MakeClass("TMiniDSTWriter"); 


TClasTool *mte2= new TClasTool();        mHtml->MakeClass("TClasTool"); 
TCounter  *count=new TCounter();         mHtml->MakeClass("TCounter"); 
TCounterCol  *cc=new TCounterCol();      mHtml->MakeClass("TCounterCol"); 
TInterrupt   *interr=new TInterrupt();   mHtml->MakeClass("TInterrupt"); 

TPartSieve *ps= new TPartSieve();       mHtml->MakeClass("TPartSieve"); 
TPartSieveHists *psh= new TPartSieveHists();       mHtml->MakeClass("TPartSieveHists"); 
TFilter *fil=new TFilter();             mHtml->MakeClass("TFilter"); 

TMapUtils *mutil=new TMapUtils();      mHtml->MakeClass("TMapUtils");
TDBItem   *dbi=new TDBItem();          mHtml->MakeClass("TDBItem"); 
TQuery    *query=new TQuery();        mHtml->MakeClass("TQuery");
TRunIndex *runi=new TRunIndex();      mHtml->MakeClass("TRunIndex");

if(gROOT->GetClass("THbookFile")){
  TNT10Reader *nt10=new TNT10Reader(); mHtml->MakeClass("TNT10Reader");
 }

 mHtml->Convert("../Utils/Write_ROOT_Tree.cc","Write_ROOT_Tree Main routine for writing DSTs","../html");

//if(!gSystem->AccessPathName("../Gsim/HTML_Load.C",kReadPermission)){
//  gROOT->Macro("../Gsim/HTML_Load.C");
//}
//
// Make an index of ALL classes. (This includes all the root classes, we 
// unfortunately do not have access to the list of classes indexed.)
//
// Note: The THtml Class does not build an index correctly if:
//   A) There is more than one class defined per .h and .cxx file.
//   B) The name of the class or class definition files contains an _ (underscore.)
//
// You thus want to avoid TMy_Favorite_Class.h and use TMyFavoriteClass instead.
//
 
cout << "\n\nMaking the Index.\n\n";

mHtml->MakeIndex();
//
// Now remove the unwanted indexes:
//
gSystem->Unlink("../html/BASE_Index.html");
gSystem->Unlink("../html/CONT_Index.html");
gSystem->Unlink("../html/ClassIndex.html");
gSystem->Unlink("../html/G3D_Index.html");
gSystem->Unlink("../html/GPAD_Index.html");
gSystem->Unlink("../html/GRAF_Index.html");
gSystem->Unlink("../html/GUI_Index.html");
gSystem->Unlink("../html/HIST_Index.html");
gSystem->Unlink("../html/HTML_Index.html");
gSystem->Unlink("../html/MATRIX_Index.html");
gSystem->Unlink("../html/META_Index.html");
gSystem->Unlink("../html/MINUIT_Index.html");
gSystem->Unlink("../html/MYSQL_Index.html");
gSystem->Unlink("../html/NET_Index.html");
gSystem->Unlink("../html/PHYSICS_Index.html");
gSystem->Unlink("../html/POSTSCRIPT_Index.html");
gSystem->Unlink("../html/RINT_Index.html");
gSystem->Unlink("../html/TREE_Index.html");
gSystem->Unlink("../html/UNIX_Index.html");
gSystem->Unlink("../html/X11_Index.html");
gSystem->Unlink("../html/X11TTF_Index.html");
gSystem->Unlink("../html/METAUTILS_Index.html");
gSystem->Unlink("../html/TREEPLAYER_Index.html");
gSystem->Unlink("../html/HBOOK_Index.html");
//
// Rename the USER__Index.html to ClasTool_index.html, replace USER with ClasTool and
// also add the ExtraIndex.html contend right at the end, by <!--SIGNATURE-->
 if(gROOT->GetVersionInt() > 51200 ){
   gSystem->Exec("sed 's/USER_/ClasTool/;/<!--SIGNATURE-->/rExtraIndex' ../html/USER__Index.html > ../html/ClasTool_index.html");
   gSystem->Unlink("../html/USER__Index.html");
 }else{
   gSystem->Exec("sed 's/USER/ClasTool/;/<!--SIGNATURE-->/rExtraIndex' ../html/USER_Index.html > ../html/ClasTool_index.html");
   gSystem->Unlink("../html/USER_Index.html");
 }

gSystem->Exit(0); // Exit without error, so make will be happy.
}


