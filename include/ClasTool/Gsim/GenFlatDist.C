//
// This is really the dumbest thing you can do, and because of that, not all that bad.
//
// Produce an electron, and a proton pretty much randomly.
// The data is thrown in between limits for omega and theta_e, and final proton momentum and theta_p.
// 
// If you wish to *COMPILE* this code, and run it as a loaded shared library follow the following steps:
// 1) Uncomment the following line:
// #include "TGsim.h"
// 2) Load needed libs, at root prompt type:
//    gSystem->Load("libE2Banks.so");gSystem->Load("libE2AnaTool");gSystem->Load("libPhysics.so");gSystem->Load("libGsim.so");
// 3) Set the include path, at root prompt type:
//    gSystem->SetIncludePath("-I$ROOTSYS/include -I/home/maurik/codes/clas/release-3-0/packages/bosio/ -I/home/maurik/e2run/root/E2AnaTool -I/home/maurik/e2run/root/e2banks");
// 4) Compile and load:
//    .L ThrowFlat.C++
//
// NOTE: The ++ at the end, that causes root to compile the code first, so it is NOT interpreted.

void GenFlatDist(Int_t Nevent=200)
{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGsim.so");
  
  Int_t Max_Events = Nevent;

  Double_t Ebeam = 2.2;
  Double_t Omega_min = 0.1;
  Double_t Omega_max = 2.0;
  Double_t Theta_min = 10.;
  Double_t Theta_max = 90.;
  
  Double_t Pfinal_min=0.15;
  Double_t Pfinal_max=2.0;
  Double_t Ptheta_min=5.;
  Double_t Ptheta_max=150.;

  Double_t Vz_min=-4.;   // Min Z of target.
  Double_t Vz_max=0.25;  // Max Z of target.
  Double_t Vradius=0.25;
  TVector3 V3Voffset(0,0,0);

  TGsim *gs = new TGsim();
  gs->OpenFile("flatdist.bos","w");

  TRandom rnd;

  Double_t Vz;
  Double_t Vr;
  Double_t Vph;
  Double_t Omega;
  Double_t Pe;
  Double_t Theta_e;
  Double_t Phi_e;
  Double_t Pf;
  Double_t Theta_p;
  Double_t Phi_p;

  TVector3 V3vert(0,0,0);
  TLorentzVector V4mom(0,0,0,0);

  for(Int_t i=0;i<Max_Events;i++){

    gs->FillHEAD(1,i);

    //    gs->CreateBank("PART",0,-1,2); // Space for 2 particles.
    gs->CreateBank("MCTK",0,11,2); // Space for 2 particles.
    gs->CreateBank("MCVX",0,5,2); // Space for 2 particles.

    TVector3  v3_prot,v3_pim;

    double  phi_angle_prot = rnd.Rndm()*TMath::Pi()*2.;

    v3_prot.SetMagThetaPhi(1.5,50./57.29,phi_angle_prot);
    v3_pim.SetMagThetaPhi(1.7,65./57.29,TMath::Pi()-phi_angle_prot);
    
    TVector3  part_vrtx(0.,0.,-10.);

    V4mom.SetVectM(v3_prot,0.938);

    gs->FillMCTKPART(0,2212,1,part_vrtx,V4mom);
    V4mom.SetVectM(v3_pim,0.135);
    gs->FillMCTKPART(1,-211,-1,part_vrtx,V4mom);



    //    gs->PrintPart();

    gs->WriteBOS();
  }

  cerr << "Calling file close \n.";

  gs->CloseFiles();

  cerr << "File close done\n.";

}
