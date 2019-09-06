{
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasBanks.so");
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasTool.so");
TMiniDSTReader *fReader = new TMiniDSTReader();

fReader->AddFile("/scratch/clas/e1d/schemed_root_old/schemed_rdst23456_13.root");
Int_t iCount  = 0;
Int_t iECount = 0;
while( fReader->Next()==0){
  //   fReader->PrintEventStatus();
  if(fReader->GetNPart()>0) iECount++;
  iCount++;
   //   cout << fReader->GetNRows() << endl;
}
cout << " Count = " << iCount << " G evt " << iECount << endl;
}
