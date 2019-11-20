#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "target_Info.h"
#include "e2a_constants.h"
#include "constants.h"
//#include "Fiducial.h"


// NOTE: Not doing anything with Vtx Cor YET
// NOTE: You can input any nucleus to scatter off of, but it should match the ones the SRC and Mean Field generator was run off of
// The purpose of this file is to add smearing, acceptances, and fiducial cuts to pseudodata for the Clas12 E2A reactor

int main(int argc, char** argv){
  if (argc != 5){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try ./simulator [Mean Field (1) or SRC (2)] [Atomic Number] [input = ../../SRC_Project/Build/DATA/gen_file.root] [output = Sam_analysis/DATA/simulator_(SRC or MF).root]\n";
    return -1;
  }

// TTree to save values during the loop
  TFile * input_file = new TFile(argv[3]);          // MF or SRC input file
  TFile * outfile = new TFile(argv[4],"RECREATE");  // New file with smearing, acceptances, and fiducials
  TTree * generator_tree;                           // MF or SRC pseudodata tree 
  TTree * outtree;                                  // New tree with smearing, acceptances, and fiducials
  // take in the correct input/output tree for Mean Field data
  if (atoi(argv[1]) == 1.){
    std::cout << "Performing smearing, acceptances, and fiducial cuts for Mean Field Data \n";
    generator_tree = (TTree*)input_file->Get("genT");
    outtree = new TTree("T", "Mass Tree");}
  // take in the correct input/output tree for SRC data
  else if (atoi(argv[1]) == 2.){
    std::cout << "Performing smearing, acceptances, and fiducial cuts for SRC Data \n";
    generator_tree = (TTree*)input_file->Get("genT");
    outtree = new TTree("T", "SRC Tree");}
  // If none specified
  else{
    std::cout << "Input 1 for Mean Field data and 2 for SRC data\n";
    return -1;}
    

// Define New Tree Variables for massT
  const int maxPart = 50.;  // arbitrary length to store each ejected particle in array
  double mom_x[maxPart], mom_y[maxPart], mom_z[maxPart], vtx_z_cor[maxPart];
  int nParticles, Part_type[maxPart];
  double Xb, weighted, Q2;
  TVector3 q_vec;
  
//Variables we are saving (match skim_tree in E2A): "name", &variable, "name with elements"/double
  outtree->Branch("nParticles", &nParticles, "nParticles/I"            ); // number of particles detected (ejected), integer -> 2 for MF, 3 for SRC  
  outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" ); // neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" ); // vertext Z corrected, arrays of integer, length of nParticles
  outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     ); // array of final momentum of ejected particles in x direction, length of nParticles, GeV
  outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     ); // array of final momentum of ejected particles in z direction, length of nParticles, GeV
  outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     ); // array of final momentum of ejected particles in y direction, length of nParticles, GeV
  outtree->Branch("Xb"        , &Xb        , "Xb/D"                    ); // Bjorken X
  outtree->Branch("Q2"        , &Q2        , "Q2/D"                    ); // Q^2
  outtree->Branch("weighted"  , &weighted , "weighted/D"               ); // weighted differential cross section
  

// Define Old Tree Variables I need from genT
  double QSq, weight, p_acc, e_acc;
  double pLead[3],  pRec[3], pe[3];
  int lead_type, rec_type;  
  
// Get branch (values) from MF_generator.cpp (MEAN FIELD) or gen_weight.cpp (SRC)
  // Mean Field trees
  if (atoi(argv[1]) == 1.){
    generator_tree->SetBranchAddress("lead_type",&lead_type);  // type of ejected nucleon
    generator_tree->SetBranchAddress("pLead",pLead);           // momentum vector of ejected nucleon
    generator_tree->SetBranchAddress("pe",pe);                 // momentum vector of scattered electron (final)
    generator_tree->SetBranchAddress("weight",&weight);        // weighted differential cross section
    pRec[0] = pRec[1] = pRec[2] = 0.001;}                      // deleted at end of code-> Placeholder for math
  // SRC trees
  else if (atoi(argv[1]) == 2.){
    generator_tree->SetBranchAddress("lead_type",&lead_type);  // type of ejected nucleon (1rst nucleon in SRC pair)
    generator_tree->SetBranchAddress("rec_type",&rec_type);    // type of ejected nucleon (2nd nucleon in SRC pair)
    generator_tree->SetBranchAddress("pLead",pLead);           // momentum vector of ejected nucleon (1rst nucleon in SRC pair)
    generator_tree->SetBranchAddress("pRec",pRec);             // momentum vector of ejected nucleon (2nd nucleon in SRC pair)
    generator_tree->SetBranchAddress("pe",pe);                 // momentum vector of scattered electron (final)
    generator_tree->SetBranchAddress("weight",&weight);}       // weighted differential cross section

  
//Initialize CLasses:
    const int A = atoi(argv[2]);  // Takes in the atomic number of the nucleus we scattered off of
    TRandom3 myRand(0);           // Initializes a random number generator
    target_Info find_acc(A);      // Initializes the acceptance and fiducial class in target_info.cpp
    // Fiducial *targFid = new Fiducial(4461,2250,5996,"4He", true);  // Beam energy, Magnetic Field Strength, Minature Magnetic Field Strength, Target_String, ??
    
// Loop over all entries
  for (int i = 0; i < generator_tree->GetEntries(); i++){   
    generator_tree->GetEvent(i);
    vtx_z_cor[1] = -1.;  // tbd later on
    p_acc = e_acc = 0;   // Set to zero after every loop
    TVector3 pe_TVec(pe[0], pe[1], pe[2]);              // TVector3 of the final electron momentum 
    TVector3 pLead_TVec(pLead[0], pLead[1], pLead[2]);  // TVector2 of final nucleon momentum (1rst nucleon in SRC pair)
    TVector3 pRec_TVec(pRec[0], pRec[1], pRec[2]);      // TVector2 of final nucleon momentum (2nd nucleon in SRC pair)
    
    // Assign variables: nucleon type, number of particles ejected
    Part_type[0] = eCode;      // particle code for electron
    Part_type[1] = lead_type;  // particle code for scattered nucleon (proton or neutron)
    if (atoi(argv[1]) == 1.){
      nParticles = 2.;}        //specific for Mean Field Generator - scattered electron and nucleon}
    else if (atoi(argv[1]) == 2.){
      Part_type[2] = rec_type;      // scattered electron, scattered nucleon (proton or neutron), scattered pair nucleon
      nParticles = 3.;}              //specific for SRC Generator - scattered electron and lead + recoil nucleon}

    
    // Acceptance Cut
    if (Part_type[0] = -11.){
      e_acc = find_acc.e_acc(pe_TVec);}
    if (Part_type[1] == 2212.){ // for now we are only looking at the lead nucleon (semi-inclusive)
      p_acc = find_acc.p_acc(pLead_TVec);}
    else if (Part_type[1] == 2112.){
      p_acc = 0;}    // for now we are only doing proton ejection, need to change later if want neutron
    if (e_acc == 0. or p_acc == 0.) continue;


    // Smeared Variables (changing) 
    TVector3 pe_TVec_smear(pe[0], pe[1], pe[2]);
    pe_TVec_smear.SetMag(pe_TVec.Mag() * (1 + myRand.Gaus(0, 0.02)) );
    TVector3 pLead_TVec_smear(pLead[0], pLead[1], pLead[2]);
    pLead_TVec_smear.SetMag(pLead_TVec.Mag() * (1 + myRand.Gaus(0, 0.02)) );
    // for both MF and SRC -> really just a dummy variable if MF
    TVector3 pRec_TVec_smear(pRec[0], pRec[1], pRec[2]);
    pRec_TVec_smear.SetMag(pRec_TVec.Mag() * (1 + myRand.Gaus(0, 0.0125)) );
    

    // electrons
    mom_x[0] = pe_TVec_smear.X();
    mom_y[0] = pe_TVec_smear.Y();
    mom_z[0] = pe_TVec_smear.Z();
    // Lead
    mom_x[1] = pLead_TVec_smear.X();
    mom_y[1] = pLead_TVec_smear.Y();
    mom_z[1] = pLead_TVec_smear.Z();
    // Now only if SRC Pair do we record the second recoil nucleon
    if (atoi(argv[1]) == 2.){
      mom_x[2] = pRec_TVec_smear.X();
      mom_y[2] = pRec_TVec_smear.Y();
      mom_z[2] = pRec_TVec_smear.Z();}
    // Math to find Q2 and Xb
    q_vec[0] = -mom_x[0];            // X-momentum transfered to nucleon from photon; Gev
    q_vec[1] = -mom_y[0];            // Y-momentum transfered to nucleon from photon; Gev
    q_vec[2] = vBeam.Z() - mom_z[0]; // Z-momentum transfered to nucleon from photon; Gev
    double Ek = sqrt(me*me + pe_TVec.Mag2()); // energy of outgoing electron
    double w = vBeam.Z() - Ek;                // Energy of scattered electron; Gev
    Q2 = q_vec.Mag2() - w*w;
    if (Part_type[1] == 2212.){ // only for lead
      Xb = Q2/(2*mP*w);}
    else if (Part_type[1] == 2112.){ // only for lead
      Xb = Q2/(2*mN*w);}

    
    // Fiducial Cut    
    //if (not find_acc.pass_semi_fid(pe_TVec_smear, pLead_TVec_smear)) continue; // fiducial cut of smears electron, lead nucleon

    
    // Find New Weights
    if (p_acc > 0.75 and e_acc > 0){
      weighted = weight/(p_acc*e_acc);}
    else{
      weighted = 0;}
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







// Taking out for now:
    /**
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
    **/
