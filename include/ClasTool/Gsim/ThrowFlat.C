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

void ThrowFlat(Int_t Nevent=100)
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

    gs->FillHEAD(10,i);

    gs->CreateBank("PART",0,-1,2); // Space for 2 particles.

    Vz = Vz_min + rnd.Uniform(Vz_max-Vz_min);
    Vr = rnd.Rndm()*Vradius;
    Vph= rnd.Rndm()*TMath::Pi()*2.;

    V3vert.SetX(Vr*TMath::Cos(Vph) + V3Voffset.X());
    V3vert.SetY(Vr*TMath::Sin(Vph) + V3Voffset.Y());
    V3vert.SetZ(Vz + V3Voffset.Z());

    Omega   = Omega_min + rnd.Uniform(Omega_max-Omega_min);
    Theta_e = (Theta_min + rnd.Uniform(Theta_max-Theta_min))*TMath::Pi()/180.;
    Phi_e   = rnd.Rndm()*TMath::Pi()*2.;
    
    Pe=Ebeam-Omega;
    V4mom.SetT(Pe);
    V4mom.SetX(Pe*TMath::Sin(Theta_e)*TMath::Cos(Phi_e));
    V4mom.SetY(Pe*TMath::Sin(Theta_e)*TMath::Sin(Phi_e));
    V4mom.SetZ(Pe*TMath::Cos(Theta_e));

    gs->FillPART(0,0,3,-1,&V3vert,&V4mom);

    Pf  = Pfinal_min + rnd.Uniform(Pfinal_max-Pfinal_min);
    Theta_p = (Ptheta_min + rnd.Uniform(Ptheta_max-Ptheta_min))*TMath::Pi()/180.;
    Phi_p   = rnd.Rndm()*TMath::Pi()*2.;
    V4mom.SetT(Pf);
    V4mom.SetX(Pf*TMath::Sin(Theta_p)*TMath::Cos(Phi_p));
    V4mom.SetY(Pf*TMath::Sin(Theta_p)*TMath::Sin(Phi_p));
    V4mom.SetZ(Pf*TMath::Cos(Theta_p));

    gs->FillPART(0,1,14,1,&V3vert,&V4mom);

    //    gs->PrintPart();

    gs->WriteBOS();
  }

  cerr << "Calling file close \n.";

  gs->CloseFiles();

  cerr << "File close done\n.";

}
