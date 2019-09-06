// Read in a GSIM INPUT file, and make some simple plots to 
// check what is going in.

TGsim  *gsim;
TH1F   *elec_mom;
TH1F   *elec_tht;
TH1F   *elec_phi;
TH2F   *emomtht;
TH2F   *ethtphi;
TH1F   *prot_mom;
TH1F   *prot_tht;
TH1F   *prot_phi;
TH2F   *pmomtht;
TH2F   *pthtphi;

void ReadGSIMIn(char *file="gsim.bos",Int_t NEvents=10000){

  Int_t nevt=0;

  gsim=new TGsim();
  gsim->OpenFile(file);

  elec_mom = new TH1F("elec_mom","Electron Momentum",100,0.,2.);
  elec_tht = new TH1F("elec_tht","Electron Theta",180,0.,180.);
  elec_phi = new TH1F("elec_phi","Electron Phi",180,0.,360.);
  emomtht  = new TH2F("emomtht","Electron Theta vs. Momentum",50,0.,2.,30,0.,60.);
  ethtphi  = new TH2F("ethtphi","Electron Theta vs. Phi",60,0.,360.,60.,0.,60.);
  prot_mom = new TH1F("prot_mom","Proton Momentum",100,0.,2.);
  prot_tht = new TH1F("prot_tht","Proton Theta",180,0.,180.);
  prot_phi = new TH1F("prot_phi","Proton Phi",180,0.,360.);
  pmomtht  = new TH2F("pmomtht","Proton Theta vs Momentum",60,0.,60.,50,0.,1.);
  pthtphi  = new TH2F("pthtphi","Proton Theta vs. Phi",60,0.,60.,60.,0.,360.);


  TMCVXClass *mcvx;
  TMCTKClass *mctk;
  void *tmpptr;
  while(!gsim->ReadBOS() && ++nevt <= NEvents){
    if( nevt%1000 == 0 ){
      cerr << "Now at " << nevt << endl;
    }

    if(gsim->IsMCEvent(-4)){
    
      if( (tmpptr = gsim->GetBank("MCTK",0)) == NULL){
	cerr << "No MCTK bank at event\n";
	gsim->PrintHead();
	mctk=NULL;
      }else{
	mctk=new TMCTKClass(tmpptr,0);	
      }

      if( (tmpptr = gsim->GetBank("MCVX",0)) == NULL){
	cerr << "No MCVX bank at event\n";
	gsim->PrintHead();
	mcvx=NULL;
      }else{
	mcvx=new TMCVXClass(tmpptr,0);
      }

      if(mctk && mcvx){

	Float_t etht=180*TMath::ACos(mctk->Cz)/TMath::Pi();
	Float_t ephi=180*TMath::ASin(mctk->Cy)/TMath::Pi();
	if(mctk->Cx<0 )ephi= 180-ephi;
	if(mctk->Cx>0 && ephi < 0)ephi= 360+ephi;
	
	cerr << "Pmom: " << mctk->Pmom << " Etht: " << etht << " EPhi: " << ephi;
	cerr << " Cz: " << mctk->Cx << " Cy: " << mctk->Cy << " Cz: " << mctk->Cz << endl;
	
	elec_mom->Fill(mctk->Pmom);
	elec_tht->Fill(etht);
	elec_phi->Fill(ephi);
	emomtht->Fill(etht,mctk->Pmom);
	ethtphi->Fill(ephi,etht);
	
      }
    }
  }

}
