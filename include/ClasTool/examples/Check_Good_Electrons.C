{
//
// Quick and dirty good electron check.
//
// 
THEADERClass *head;
TEVNTClass   *evnt;
TECPBClass   *ecpb;
TDCPBClass   *dcpb;
TCCPBClass   *ccpb;
TDSTReader   *R;

Int_t No_Electron_id_Count=0;
Int_t Z_Cut_fail = 0;
Int_t Events_left=0;
Int_t Bad_SC_stat=0;
Int_t Bad_CC_stat=0;
Int_t Bad_EC_stat=0;
Int_t CC_cut=0;
Int_t EC_cut=0;

TH1F *Z_vertex;
TH1F *CCph;
TH1F *ECtot;
TH1F *ECin;
TH1F *ECout;
TH2F *ECTI;
TH2F *ECIO;
TH1F *Emom;

TH1F *C_Z_vertex;
TH1F *C_CCph;
TH1F *C_ECtot;
TH1F *C_ECin;
TH1F *C_ECout;
TH2F *C_ECTI;
TH2F *C_ECIO;
TH1F *C_Emom;
 
if(!Init)Z_vertex=new TH1F("Z_vertex","Z Vertex",200,-10,5);
if(!Init)CCph = new TH1F("CCph","CC total",200,0,200);
if(!Init)ECtot = new TH1F("ECtot","EC total",200,0,2);
if(!Init)ECin = new TH1F("ECin","EC inner",200,0,2);
if(!Init)ECout = new TH1F("ECout","EC outer",200,0,2);
if(!Init)ECTI = new TH2F("ECTI","EC tot vs EC in",100,0,2,100,0,2);
if(!Init)ECTO = new TH2F("ECTO","EC tot vs EC out",100,0,2,100,0,2);
if(!Init)ECIO = new TH2F("ECIO","EC in vs EC out",100,0,2,100,0,2);
if(!Init)ECTP = new TH2F("ECTP","EC tot vs Momentum",100,0,5,100,0,2);
if(!Init)Emom = new TH1F("Emom","Electron Momentum",200,0,5);

if(!Init)C_Z_vertex=new TH1F("C_Z_vertex","Z Vertex",200,-10,5);
if(!Init)C_CCph = new TH1F("C_CCph","CC total",200,0,200);
if(!Init)C_ECtot = new TH1F("C_ECtot","EC total",200,0,2);
if(!Init)C_ECin = new TH1F("C_ECin","EC inner",200,0,2);
if(!Init)C_ECout = new TH1F("C_ECout","EC outer",200,0,2);
if(!Init)C_ECTI = new TH2F("C_ECTI","EC tot vs EC in",100,0,2,100,0,2);
if(!Init)C_ECTO = new TH2F("C_ECTO","EC tot vs EC out",100,0,2,100,0,2);
if(!Init)C_ECIO = new TH2F("C_ECIO","EC in vs EC out",100,0,2,100,0,2);
if(!Init)C_ECTP = new TH2F("C_ECTP","EC tot vs Momentum",100,0,5,100,0,2);
if(!Init)C_Emom = new TH1F("C_Emom","Electron Momentum",200,0,5);

if(!R)R= new TDSTReader();
if(!R->GetNFiles())R->AddFile("r32966.root");
R->Next();

//for(Int_t i=0;i<R->GetEntries();i++){
for(Int_t i=1;i<300000;i++){

  R->Next();
  head = R->fEvntHeader;

  if( (i%10000)==0){
    cout << "At " << i << " Event number " << head->NEvent << endl;
  }

  if( (i%100000)==0){
    cout << "\nDone with " << i << " events \n";
    cout << "Events wo electron id  : " << No_Electron_id_Count << "\n";
    cout << "Events failing zcut    : " << Z_Cut_fail << "\n";
    cout << "Events with bad SC stat: " << Bad_SC_stat << endl;
    cout << "Events with bad CC stat: " << Bad_CC_stat << endl;
    cout << "Events with bad EC stat: " << Bad_EC_stat << endl;
    cout << "Left with              : " << Events_left << endl << endl; 
  }


  if(R->GetNPart() == 0){
    //    cout << "No particles in event at " << i << " Event number " << head->NEvent << endl;
    //    return;
    continue;
  }

  if(R->GetPartID(0) != 11){
    //    cout << "No electron 0, at " << i << " Event number " << head->NEvent << endl;
    No_Electron_id_Count++;
    // continue;
  }

  evnt = (TEVNTClass *)R->fcaEVNTStore->At(0);

  Double_t e_mom = sqrt(evnt->Px*evnt->Px+evnt->Py*evnt->Py+evnt->Pz*evnt->Pz);

  Emom->Fill(e_mom);

  Int_t sc_stat = evnt->Scstat;
  if(sc_stat<=0){
    Bad_SC_stat++;
    continue;
  }

  Int_t cc_stat = evnt->Ccstat;
  if(cc_stat<=0){
    Bad_CC_stat++;
    continue;
  }

  Int_t ec_stat = evnt->Ecstat;
  if(ec_stat<=0){
    Bad_EC_stat++;
    continue;
  }

  ccpb = (TCCPBClass *)R->fcaCCPBStore->At(cc_stat-1);
  if(!ccpb){
    cerr << "Got CCPB problem: " << cc_stat << "\n";
    cerr << "CCPBStore " << R->fcaCCPBStore->GetEntries() << endl;
    return;
    continue;
  }

  CCph->Fill(ccpb->Nphe);

  if(ccpb->Nphe < 4.){
    CC_cut++;
    continue;
  }


  ecpb = (TECPBClass *)R->fcaECPBStore->At(ec_stat-1);
  if(!ecpb){
    cerr << "Got ECPB problem: " << ec_stat << "\n";
    cerr << "ECPBStore " << R->fcaECPBStore->GetEntries() << endl;
    return;
    continue;
  }

  Z_vertex->Fill(evnt->Z);
  if(evnt->Z < -6. || evnt->Z > -2.5 ){
    Z_Cut_fail++;
    continue;
  }

  Double_t e_mom = sqrt(evnt->Px*evnt->Px+evnt->Py*evnt->Py+evnt->Pz*evnt->Pz);

  ECtot->Fill(ecpb->Etot);
  ECin->Fill(ecpb->Ein);
  ECout->Fill(ecpb->Eout);
  ECTI->Fill(ecpb->Ein,ecpb->Etot);
  ECTO->Fill(ecpb->Eout,ecpb->Etot);
  ECIO->Fill(ecpb->Eout,ecpb->Ein);
  ECTP->Fill(e_mom,ecpb->Etot);



  if(ecpb->Etot< 0.2 || ecpb->Etot < (e_mom*0.3 - 0.12) ){
    EC_cut++;
    continue;
  }

  C_Emom->Fill(e_mom);
  C_Z_vertex->Fill(evnt->Z);
  C_CCph->Fill(ccpb->Nphe);
  C_ECtot->Fill(ecpb->Etot);
  C_ECin->Fill(ecpb->Ein);
  C_ECout->Fill(ecpb->Eout);
  C_ECTI->Fill(ecpb->Ein,ecpb->Etot);
  C_ECTO->Fill(ecpb->Eout,ecpb->Etot);
  C_ECIO->Fill(ecpb->Eout,ecpb->Ein);
  C_ECTP->Fill(e_mom,ecpb->Etot);


  Events_left++;

}

    cout << "\nDone with " << i << " events \n";
    cout << "Events wo electron id  : " << No_Electron_id_Count << "\n";
    cout << "Events failing zcut    : " << Z_Cut_fail << "\n";
    cout << "Events with bad SC stat: " << Bad_SC_stat << endl;
    cout << "Events with bad CC stat: " << Bad_CC_stat << endl;
    cout << "Events with bad EC stat: " << Bad_EC_stat << endl;
    cout << "Left with              : " << Events_left << endl; 


}
