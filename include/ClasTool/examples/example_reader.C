{
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasBanks.so");
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasTool.so");
TDSTReader *fReader = new TDSTReader();
fReader->AddFile("test2.root");

for(int j =0; j < 4; j++){
   cout << fReader->Next() << endl;; // Reads next event from file..
   fReader->PrintEventStatus();
    cout << fReader->GetNRows() << endl;
}

delete fReader;
}
