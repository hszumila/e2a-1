#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "target_Info.h"


int main(int argc, char** argv){
  if (argc != 4){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try ./pseudo_skim_tree [Mean Field (1) or SRC (2)] ../DATA/gen_file.root ../DATA/pseudo_skim_tree.root\n";
    return -1;
  }

// TTree to save values during the loop
  TFile * input_file = new TFile(argv[2]);
  TFile * outfile = new TFile(argv[3],"RECREATE");
  TTree * generator_tree;
  TTree * outtree;
  if (atoi(argv[1]) == 1.){
    std::cout << "Smearing for Mean Field Data \n";
    generator_tree = (TTree*)input_file->Get("genT");
    outtree = new TTree("T", "Mass Tree");}
  else if (atoi(argv[1]) == 2.){
    std::cout << "Smearing for SRC Data \n";
    generator_tree = (TTree*)input_file->Get("genT");
    outtree = new TTree("T", "SRC Tree");}
  else{
    std::cout << "Input 1 for Mean Field and 2 for SRC Pair \n";
    return -1;}
    

// Define New Tree Variables for massT
  const double Pbz = 4.461; // THIS SHOULD BE BEAM ENERGY, GeV! CHANGE GENERATOR IF NOT
  const int maxPart = 50;
  const int Me = 0.0005109906; // Mass of electron; GeV
  double mom_x[maxPart], mom_y[maxPart], mom_z[maxPart], vtx_z_cor[maxPart];
  int nParticles, Part_type[maxPart];
  double Xb, weighted, Q2, Mn;
  TVector3 q_vec;
  
//Variables we are saving (match skim_tree in E2A): "name", &variable, "name with elements"/double
  outtree->Branch("nParticles", &nParticles, "nParticles/I"            ); //number of particles detected, integer - for me it will always be two  
  outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" ); //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" ); //vertext Z corrected, arrays of integer, length of nParticles
  outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     ); //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     ); //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     ); //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  outtree->Branch("Xb"        , &Xb        , "Xb/D"                    ); // Bjorken X
  outtree->Branch("Q2"        , &Q2        , "Q2/D"                    ); // Bjorken X
  outtree->Branch("weighted"  , &weighted , "weighted/D"               ); // weight  
  

// Define Old Tree Variables I need from genT
  int nucleon_type_lead, nucleon_type_recoil;
  double P1_prime_vec_lead[3],  P1_prime_vec_recoil[3], P_k_vec[3];
  double X_b, Q_2, weight, p_acc, e_acc;
  
  
// Get branch (values) from generator.cpp (MEAN FIELD) or gen_weight.cpp (SRC)
  if (atoi(argv[1]) == 1.){
    generator_tree->SetBranchAddress("X_b",&X_b);
    generator_tree->SetBranchAddress("nucleon_type",&nucleon_type_lead);
    generator_tree->SetBranchAddress("P1_prime_vec",P1_prime_vec_lead);
    generator_tree->SetBranchAddress("P_k_vec",P_k_vec);
    generator_tree->SetBranchAddress("weight",&weight);
    generator_tree->SetBranchAddress("Q_2",&Q_2);
    P1_prime_vec_recoil[0] = P1_prime_vec_recoil[1] = P1_prime_vec_recoil[2] = 0;} // to delete at end -> Placeholder for now
  else if (atoi(argv[1]) == 2.){
    generator_tree->SetBranchAddress("xB",&X_b);
    generator_tree->SetBranchAddress("lead_type",&nucleon_type_lead);
    generator_tree->SetBranchAddress("rec_type",&nucleon_type_recoil);
    generator_tree->SetBranchAddress("pLead",P1_prime_vec_lead);
    generator_tree->SetBranchAddress("pRec",P1_prime_vec_recoil);
    generator_tree->SetBranchAddress("pe",P_k_vec);    
    generator_tree->SetBranchAddress("weight",&weight);
    generator_tree->SetBranchAddress("QSq",&Q_2);}

  
//Initialize CLasses:
    const int A = 4;      // Assuming We are scattering off Helium
    TRandom3 myRand(0); // Get random number
    target_Info find_acc(A); // assumes a target of helium-4
    
// Loop over all entries
  for (int i = 0; i < generator_tree->GetEntries(); i++){   
    generator_tree->GetEvent(i);
    p_acc = e_acc = 0;
    
    // Predefined for this situation; change if changed experiment
    Part_type[0] = -11.;
    if (atoi(argv[1]) == 1.){
      Part_type[1] = nucleon_type_lead; // scattered electron, scattered nucleon (proton or neutron)
      nParticles = 2.;} //specific for Mean Field Generator - scattered electron and nucleon}
    else if (atoi(argv[1]) == 2.){
      Part_type[1] = nucleon_type_lead;   // scattered electron, scattered nucleon (proton or neutron), scattered pair nucleon
      Part_type[2] = nucleon_type_recoil; // scattered electron, scattered nucleon (proton or neutron), scattered pair nucleon
      nParticles = 3.;} //specific for SRC Generator - scattered electron and lead + recoil nucleon}
    vtx_z_cor[1] = -1.; // tbd later on

    /**
    // Variable to take from old code (unsmeared)
    Xb = X_b;
    Q2 = Q_2;
    weighted = weight;
    mom_x[0] = P_k_vec[0];
    mom_y[0] = P_k_vec[1];
    mom_z[0] = P_k_vec[2];
    mom_x[1] = P1_prime_vec[0];
    mom_y[1] = P1_prime_vec[1];
    mom_z[1] = P1_prime_vec[2];
    **/

    // Smeared Variables (changing) 
    TVector3 Pk_TVec(P_k_vec[0], P_k_vec[1], P_k_vec[2]);
    TVector3 P1_prime_TVec_lead(P1_prime_vec_lead[0], P1_prime_vec_lead[1], P1_prime_vec_lead[2]);
    TVector3 P1_prime_TVec_recoil(P1_prime_vec_recoil[0], P1_prime_vec_recoil[1], P1_prime_vec_recoil[2]);
    // lead nucleon
    double P1_prime_lead_mag = P1_prime_TVec_lead.Mag()*(1. + myRand.Gaus(0, 0.02));
    double_t theta_P1_prime_lead = P1_prime_TVec_lead.Theta();
    double_t Phi_P1_prime_lead = P1_prime_TVec_lead.Phi();
    // recoil nucleon
    double P1_prime_recoil_mag = P1_prime_TVec_recoil.Mag()*(1. + myRand.Gaus(0, 0.02));
    double_t theta_P1_prime_recoil = P1_prime_TVec_recoil.Theta();
    double_t Phi_P1_prime_recoil = P1_prime_TVec_recoil.Phi();
    // final electron
    double Pk_mag = Pk_TVec.Mag()*(1. + myRand.Gaus(0, 0.0125));
    double_t theta_Pk = Pk_TVec.Theta();
    double_t Phi_Pk = Pk_TVec.Phi();
    
    mom_x[0] = Pk_mag*sin(theta_Pk)*cos(Phi_Pk);      // Scattered electron's X-momentum; Gev  
    mom_y[0] = Pk_mag*sin(theta_Pk)*sin(Phi_Pk);      // Scattered electron's Y-momentum; Gev
    mom_z[0] = Pk_mag*cos(theta_Pk);                  // Scattered electron's Z-momentum; Gev
    TVector3 Pk_TVec_smear(mom_x[0], mom_y[0], mom_z[0]);
    mom_x[1] = P1_prime_lead_mag*sin(theta_P1_prime_lead)*cos(Phi_P1_prime_lead); // Scattered nucleon's X-momentum; Gev  
    mom_y[1] = P1_prime_lead_mag*sin(theta_P1_prime_lead)*sin(Phi_P1_prime_lead); // Scattered nucleon's Y-momentum; Gev
    mom_z[1] = P1_prime_lead_mag*cos(theta_P1_prime_lead);                        // Scattered nucleon's Z-momentum; Gev
    TVector3 P1_prime_TVec_lead_smear(mom_x[1], mom_y[1], mom_z[1]);
    // Now only if SRC Pair do we record the second recoil nucleon
    if (atoi(argv[1]) == 2.){
      mom_x[2] = P1_prime_recoil_mag*sin(theta_P1_prime_recoil)*cos(Phi_P1_prime_recoil); // Scattered nucleon's X-momentum; Gev  
      mom_y[2] = P1_prime_recoil_mag*sin(theta_P1_prime_recoil)*sin(Phi_P1_prime_recoil); // Scattered nucleon's Y-momentum; Gev
      mom_z[2] = P1_prime_recoil_mag*cos(theta_P1_prime_recoil);                          // Scattered nucleon's Z-momentum; Gev
      TVector3 P1_prime_TVec_recoil_smear(mom_x[2], mom_y[2], mom_z[2]);}

    q_vec[0] = -mom_x[0];      // X-momentum transfered to nucleon from photon; Gev
    q_vec[1] = -mom_y[0];      // Y-momentum transfered to nucleon from photon; Gev
    q_vec[2] = Pbz - mom_z[0]; // Z-momentum transfered to nucleon from photon; Gev
    double Ek = sqrt(Me*Me + Pk_mag*Pk_mag); // energy of outgoing electron
    double w = Pbz - Ek;                     // Energy of scattered electron; Gev
    

    // Get Acceptances
    if (Part_type[0] = -11){
      e_acc = find_acc.e_acc(Pk_TVec_smear);}
    for (int nucleon_ejected = 1; nucleon_ejected < atoi(argv[1])+1.; nucleon_ejected++){
      if (Part_type[1] == 2212.){
	Mn = 0.93827231; // Mass of exiting proton; Gev
	p_acc = find_acc.p_acc(P1_prime_TVec_lead_smear);}
      else if (Part_type[1] == 2112.){
	Mn = 0.93957;    // Mass of exiting neutron; Gev
	p_acc = 0;}}    // for now we are only doing proton ejection, need to change later if want neutron
    
    Q2 = q_vec.Mag2() - w*w;
    Xb = Q2/(2*Mn*w);

    //  std::cout << "p_acc: " << p_acc << "weighted1: " << weighted << "\n";
    
    // Find Acceptances
    if (p_acc != 0 and e_acc != 0){
      weighted = weight/(p_acc*e_acc);}
    else{
      weighted = 0;}

    //  std::cout << "p_acc: " << p_acc << "weighted2: " << weighted << "\n";

    // Fill Tree
      if (weighted != 0){
	outtree->Fill();}
  }


// Clean up
  input_file->Close();
  outfile->cd();
  outtree->Write();
  outfile->Close();
  return 0;
  
}






    
// Stuf from Axel <3
  /* 
     =========================
     Part_type values	(http://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html)
     -11  = electron
     2212 = proton
     2112 = neutron
     +211 = pi+
     -211 = pi-
     +111 = pi0
     =========================
  */


/** Axel's Definitions:
  outtree->Branch("nRun"      , &nRun      , "nRun/I"                  );
  outtree->Branch("nParticles", &nParticles, "nParticles/I"            );
  outtree->Branch("nProtons"  , &nProtons  , "nProtons/I"              );
  outtree->Branch("nNeutrons" , &nNeutrons , "nNeutrons/I"             );
  outtree->Branch("nPiplus"   , &nPiplus   , "nPiplus/I"               );
  outtree->Branch("nPiminus"  , &nPiminus  , "nPiminus/I"              );
  outtree->Branch("t0"        , &t0        , "t0/D"                    );
  outtree->Branch("Nu"        , &Nu        , "Nu/D"                    );
  outtree->Branch("Q2"        , &Q2        , "Q2/D"                    );
  outtree->Branch("Xb"        , &Xb        , "Xb/D"                    );
  outtree->Branch("charge"    ,  charge    , "charge[nParticles]/D"    );
  outtree->Branch("beta"      ,  beta      , "beta[nParticles]/D"      );
  outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" );
  outtree->Branch("vtx_z_unc" ,  vtx_z_unc , "vtx_z_unc[nParticles]/D" );
  outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" );
  outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     );
  outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     );
  outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     );
  outtree->Branch("phi"     ,  phi     , "phi[nParticles]/D"     );
  outtree->Branch("theta"     ,  theta     , "theta[nParticles]/D"     );
  outtree->Branch("e_deltat"  ,  e_deltat  , "e_deltat[nParticles]/D"  );
  outtree->Branch("stat_sc"   ,  stat_sc   , "stat_sc[nParticles]/I"   );
  outtree->Branch("stat_dc"   ,  stat_dc   , "stat_dc[nParticles]/I"   );
  outtree->Branch("stat_ec"   ,  stat_ec   , "stat_ec[nParticles]/I"   );
  outtree->Branch("sc_time"   ,  sc_time   , "sc_time[nParticles]/D"   );
  outtree->Branch("sc_path"   ,  sc_path   , "sc_path[nParticles]/D"   );
  outtree->Branch("sc_pad"    ,  sc_pad    , "sc_pad[nParticles]/I"    );
  outtree->Branch("ec_time"   ,  ec_time   , "ec_time[nParticles]/D"   );
  outtree->Branch("ec_path"   ,  ec_path   , "ec_path[nParticles]/D"   );
  outtree->Branch("ec_in"     ,  ec_in     , "ec_in[nParticles]/D"     );
  outtree->Branch("ec_out"    ,  ec_out    , "ec_out[nParticles]/D"    );
  outtree->Branch("ec_tot"    ,  ec_tot    , "ec_tot[nParticles]/D"    );
  outtree->Branch("ec_x"      ,  ec_x      , "ec_x[nParticles]/D"      );
  outtree->Branch("ec_y"      ,  ec_y      , "ec_y[nParticles]/D"      );
  outtree->Branch("ec_z"      ,  ec_z      , "ec_z[nParticles]/D"      );
  outtree->Branch("ec_u"      ,  ec_u      , "ec_u[nParticles]/D"      );
  outtree->Branch("ec_v"      ,  ec_v      , "ec_v[nParticles]/D"      );
  outtree->Branch("ec_w"      ,  ec_w      , "ec_w[nParticles]/D"      );
  outtree->Branch("Mass"      ,  Mass      , "Mass[nParticles]/D"      );

**/
