Int_t Load_A_Lib(char *libname){
  Int_t retval;
  
  //  cerr << "Attempting to load " << libname << endl;
  
  if(gSystem->DynamicPathName(libname,0)){
    retval=gSystem->Load(gSystem->DynamicPathName(libname,0));
    cerr << gSystem->DynamicPathName(libname,0) << " loaded. \n";
    return(0);
  }else{
    cerr << "I could not find " << libname << " in the LD_LIBRARY path,\n";
    cerr << "please check the nuber and dial again.\n";
    return(1);
  }  
}

Int_t Load_Select(){
  Int_t retval;
  retval=0;

  if(!gSystem->AccessPathName("./libSelect.so",kReadPermission)){
    retval=gSystem->Load("./libSelect.so");
    if(!retval){
      cerr << "Selections Library libSelect.so loaded from local.\n";
      cerr << SelectIdentify();
      return(0);
    }else{
      cerr << "Load_Select():: ERROR -- Could not load local libSelect.so \n";
    }
  }

  retval = Load_A_Lib("libSelect.so");
  if(retval == 0){
      cerr << "Selections Library libSelect.so loaded from LD_LIBRARY_PATH\n";
      cerr << SelectIdentify();
  }    

  return(retval);
}
