{

gROOT->Reset();

TMapTools T("quark.unh.edu","calib","dbuser");

Int_t nResult;

TFloatBuffer Energy;
T->SetDbVersion("PROD");
T->GetMapFloatBuff("RUN_CONTROL","beam","energy",18350,&Energy);
cout << "Energy: " << Energy.Get(0) << endl;

//T->SetDebugMode(0x20);
TIntBuffer Dr;
T->GetMapIntBuff("SCALERS","TRGS","values",17874,&Dr,"suleiman");

cout << " Result is array of dimmension " <<  nResult << endl;

Energy->Print();

}
