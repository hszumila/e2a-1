{
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasBanks.so");
gSystem->Load("/scratch/software/codes/slib/LinuxRH7/libClasTool.so");


TClasTool *fReader = new TClasTool();
fReader->InitDSTReader("ROOTDSTR");

fReader->AddFile("rdst23138_02.root");

Int_t head_count = 0;
Int_t evnt_count = 0;
Int_t ecpb_count = 0;
Int_t scpb_count = 0;
Int_t dcpb_count = 0;
Int_t ccpb_count = 0;
Int_t lcpb_count = 0;
Int_t tgbi_count = 0;

TH1F *hist_helicity = new TH1F("hist_helicity","helicity",3,-1.5,1.5);
while(fReader->Next() == 0){
  //   cout << fReader->Next() << endl;; // Reads next event from file..
    Int_t helicity  = fReader->GetHelicity();
    Int_t real_hel;
    Int_t id = fReader->GetPartID(0);
    int test = helicity&0x8000;
    //    cout << "ID = " << id << endl;
    if(((int)(helicity/32768)) != 0){
      real_hel = -1;
      //      cout << "Helicity = -1" << endl;
    } else {
      real_hel = 1;
      //      cout << "Helicity = 1" << endl;
    }
    cout << helicity << "   " <<  real_hel  << "  " << test << endl;
    if(id==11)
      hist_helicity->Fill(real_hel);
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
