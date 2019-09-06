//--------------------------------------------------------------------------
// This macro loads all the libraries you need to run the 
// CLASAnaTool package plus all the other packages in the
// e2run/root directory.
//
// Load this file into an interactive root session using: 
// .x rootlogin.C
//
// Alternatively you can also use gROOT->Macro("../macros/Load_Libs");
// which will run the script and load these libraries. This is useful
// in a .rootlogon.C 
//
// This script assumes that you have set your LD_LIBRARY_PATH to contain
// the libraries you are loading. A typical LD_LIBRARY_PATH can be set
// as follows:
//
// setenv LD_LIBRARY_PATH  $ROOTSYS/lib:~/slib/Linux
//
// Also, this script will work from any directory IF you set the following:
//
// setenv CLASTOOL ~/ClasTool        (replace with appropriate directory.)
//--------------------------------------------------------------------------
// Author:  Maurik Holtrop (UNH) 28/4/2000 
//
// Modified by G.Gavalian (UNH) 04/01/2002 
// to be used with new TCLASTOOL Package
//

{
  Int_t retval;

  TString gEnvCLASTOOL;
  TString gEnvLibPath;
  TString Test;

  gEnvCLASTOOL= gSystem->Getenv("CLASTOOL");
  gEnvLibPath= gSystem->Getenv("LD_LIBRARY_PATH");

  // CINT does not allow a funcion to be defined in the executing body.
  // It must be loaded, so load functions here.

  if(gEnvCLASTOOL.Length()){
    Test=gEnvCLASTOOL;
    Test.Append("/macros/Helper_Load_Libs.h");
    if(!gSystem->AccessPathName(Test,kReadPermission)){
      gROOT->LoadMacro(Test);
    }else{
      cerr << "CLASTOOL environment variable set but macros not found in: \n";
      cerr << gEnvROOT << endl;
      gSystem->Exit(1);
    }
  }else{    
    if(!gSystem->AccessPathName("../macros/Helper_Load_Libs.h",kReadPermission)){
      gROOT->LoadMacro("../macros/Helper_Load_Libs.h");
    }else if(!gSystem->AccessPathName("macros/Helper_Load_Libs.h",kReadPermission)){
      gROOT->LoadMacro("macros/Helper_Load_Libs.h");
    }else if(!gSystem->AccessPathName("Helper_Load_Libs.h",kReadPermission)){
      gROOT->LoadMacro("Helper_Load_Libs.h");
    }else{
      cerr << "Could not load helper functions from Load_Libs_Helper.h -- ABORT \n" ;
      cerr << "Please set CLASTOOL environment variable, i.e.: setenv CLASTOOL ~/ClasTool \n";
      gSystem->Exit(1);
    }
  }

  cerr<<"Welcome to the CLAS Root Analysis framework\n\n";

  // NOTE: We are now relying on dynamic loading of all needed libraries by
  // including them as requirements in the libraries themselves. See Makefiles.

  // Load the libPhysics.so which is not yet a standard in ROOT.
  //
  // if(gSystem->DynamicPathName("libPhysics.so",0) == 0){
  // cerr << "Could not find the libPhysics.so library, which should be part or ROOT.\n" ;
  // cerr << "Please check your ROOT installation, and make sure LD_LIBRARY_PATH is\n";
  // cerr << "set correctly.\n\n";
  // exit(1);
  // }else{
  // if(!(retval=gSystem->Load(gSystem->DynamicPathName("libPhysics.so",0))))
  // cerr << "Library libPhysics.so loaded: " << retval << endl;
  // else{
  // cerr << "Error loading libPhysics.so or already loaded: " << retval << endl;
  //      exit(1);
  //    }
  //  }

  Test=gEnvCLASTOOL;
  Test.Append("/include/TClasTool_Const.h");    
  if (gEnvCLASTOOL.Length() &&
      !gSystem->AccessPathName(Test,kReadPermission)){
    gROOT->LoadMacro(Test);
    cout << "Found the constants from TClasTool_Const.h" << endl;
  }else if ( !gSystem->AccessPathName("../TClasTool/TClasTool_Const.h",kReadPermission)){
    gROOT->LoadMacro("../CLASTool/TCLASTool_Const.h");
    cout << "Found the constants from TCLASTool_Const.h in TClasTool" << endl;
  }

  //  if(Load_A_Lib("libClasBanks.so")){exit(1);}
  if(Load_A_Lib("libClasTool.so")){exit(1);}

  if(gSystem->Getenv("MYSQL_INCLUDE")){
    if(Load_A_Lib("libMapUtils.so")){cerr << "No MapUtils compiled.";}
  }

}


