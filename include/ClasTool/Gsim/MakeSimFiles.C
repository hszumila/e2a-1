
void MakeSimFiles(Int_t Nevent=200)
{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGsim.so");
  gSystem->Load("libVirtualReader.so");
  gSystem->Load("libCTEventSelector.so");

  char infname[128];
  char dtline[4096];

  TCTTextReader fTR;

  TGsim *gs = new TGsim();
  gs->OpenFile("flatdist.bos","w");

  fTR.Open("gen_phi.txt");
  //  ifstream  idata(infname);
  int  MaxEvt = 50000;

  for(Int_t i=0;i<MaxEvt;i++){

    //idata.getline(dtline,4096);
    if(i%20000==0) cout << i << " events processed " << endl;
    fTR.Next();
    

    gs->FillHEAD(1,i);

    //    gs->CreateBank("PART",0,-1,2); // Space for 2 particles.

    
    gs->CreateBank("MCTK",0,11,4); // Space for 2 particles.
    gs->CreateBank("MCVX",0,5,4); // Space for 2 particles.

    TVector3  v3_prot,v3_pim,v3_pip;
    TVector3  part_vrtx, part_vect;
    TLorentzVector V4mom;

    fTR->GetNumber(0);
    part_vrtx.SetXYZ(0.,0.,0.);
    part_vect.SetXYZ(0.,0.,0.);
    V4mom.SetVectM(part_vect,0.);
    gs->FillMCTKPART(0,22,1,part_vrtx,V4mom);
    par_vect = fTR.GetVect(1);
    par_vrtx = fTR.GetVect(4);
    V4mom.SetVectM(part_vect,0.938);
    gs->FillMCTKPART(1,2212,1,part_vrtx,V4mom);
    par_vect = fTR.GetVect(7);
    par_vrtx = fTR.GetVect(10);
    V4mom.SetVectM(part_vect,0.1396);
    gs->FillMCTKPART(2,211,1,part_vrtx,V4mom);
    par_vect = fTR.GetVect(13);
    par_vrtx = fTR.GetVect(16);
    V4mom.SetVectM(part_vect,0.1396);
    gs->FillMCTKPART(3,-211,-1,part_vrtx,V4mom);

    //    gs->PrintPart();
    gs->WriteBOS();
  }

  cerr << "Calling file close \n.";

  gs->CloseFiles();

  cerr << "File close done\n.";

}
