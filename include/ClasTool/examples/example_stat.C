{
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasBanks.so");
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasTool.so");
TDSTReader *fReader = new TDSTReader();
fReader->AddFile("rdst23138_02.root");

Int_t head_count = 0;
Int_t evnt_count = 0;
Int_t ecpb_count = 0;
Int_t scpb_count = 0;
Int_t dcpb_count = 0;
Int_t ccpb_count = 0;
Int_t lcpb_count = 0;
Int_t tgbi_count = 0;


while(fReader->Next() == 0){
  //   cout << fReader->Next() << endl;; // Reads next event from file..
  //   fReader->PrintEventStatus();
   evnt_count += fReader->GetNRows("EVNT");
   ecpb_count += fReader->GetNRows("ECPB");
   scpb_count += fReader->GetNRows("SCPB");
   dcpb_count += fReader->GetNRows("DCPB");
   ccpb_count += fReader->GetNRows("CCPB");
   lcpb_count += fReader->GetNRows("LCPB");
   tgbi_count += fReader->GetNRows("TGBI");

   //    cout << fReader->GetNRows() << endl;
    
}

cout << " EVNT " << evnt_count << endl;
cout << " ECPB " << ecpb_count << endl;
cout << " SCPB " << scpb_count << endl;
cout << " DCPB " << dcpb_count << endl;
cout << " CCPB " << ccpb_count << endl;
cout << " LCPB " << lcpb_count << endl;
cout << " TGBI " << tgbi_count << endl;

delete fReader;
}
