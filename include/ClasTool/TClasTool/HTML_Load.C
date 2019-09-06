{
SourceDirs.Append(":../TClasTool"); // Tell root where the sources are
mHtml->SetSourceDir(SourceDirs);  
// Load the lib.
if(Load_A_Lib("libClasTool.so")){cerr << "ERROR loading.\n";gSystem->Exit(1);}
TClasTool *mte2= new TClasTool();        mHtml->MakeClass("TClasTool"); 
TCounter  *count=new TCounter();         mHtml->MakeClass("TCounter"); 
TCounterCol  *cc=new TCounterCol();      mHtml->MakeClass("TCounterCol"); 
TInterrupt   *interr=new TInterrupt();   mHtml->MakeClass("TInterrupt"); 
//
// Now make an html copy of the examples also.
//
//mHtml->Convert("example1.C","Simple Example how to use the TE2AnaTool Classes");
//mHtml->Convert("example2.C","More complete Example how to use the TE2AnaTool Classes");
//mHtml->Convert("TE2AnaTool_Const.h","Constants and Enumerations used in TE2AnaTool.");
}
