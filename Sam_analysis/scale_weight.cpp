#include <iostream>
#include <cmath>

#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"


// The purpose of this code is to scale the SRC and Mean Field DATA to match the Real DATA.
// Pick weighting range as it will be applied to the total data set
// The reaosn for not finding the weight over the whole set is that some data has poor events in certain regions (bias)
// We are also plotting the P_Miss distribution (as it is a well known experimental result)

int main(int argc, char** argv){
  bool applied_simulator;
  if (argc > 6. or argc < 5.){
    std::cout << "Wrong number of arguments\n";
    std::cout << "Try ./scale_weight [input MF Root DATA] [input SRC Root Data] [input Real Root Data (e2a_4GeV_He4_all.root)] [outfile new weights for SRC + MF] [Applied Simulator.cpp already: type anything; else: skip]\n";
    std::cout << "NOTE: MF and SRC input cannot be the same file.\n";
    return -1;
  }
  if (argc == 5){
    std::cout << "Using variable names assuming no simulator\n";
    applied_simulator = false;}
  if (argc == 6.){
    std::cout << "Using variable names assuming simulator was applied\n";
    applied_simulator = true;}


// Variables to Change in Program
  const int total_bins = 20.;           // Total bins for P_miss histograms
  const int lower_bound = 50.;          // P_miss axis
  const int upper_bound = 600.;         //P_miss axis
  const int start_p_scale = 300.;       // Initial P_miss value to consider in scaling (DATA comparison between real data)
  const int end_p_scale = upper_bound;  // Final P_miss calue to consider in scaling (DATA comparison between real data)
  const double xB_max_cut = 2.;         // xB cut for data
  const double xB_min_cut = 1.;         // xB cut for data
  const double QSq_max_cut = 4.1;       // Q2 cut for data
  const double QSq_min_cut = 0.5;       // Q2 cut for data
  
// Make Histograms of Missing Momentum to Compare
  TH1D * P_miss_REAL_hist = new TH1D("P_miss From REAL DATA","Missing Momentum [MeV]; Counts", total_bins,  lower_bound, upper_bound);
  TH1D * P_miss_SRC_hist = new TH1D("P_miss From SRC DATA","Missing Momentum [MeV]; Counts", total_bins,  lower_bound, upper_bound);
  TH1D * P_miss_MF_hist = new TH1D("P_miss From MF DATA","Missing Momentum [MeV]; Counts", total_bins,  lower_bound, upper_bound);
  P_miss_REAL_hist->Sumw2();
  P_miss_SRC_hist->Sumw2();
  P_miss_MF_hist->Sumw2();

// Get Information from Real DATA
  TFile * input_file_Real_DATA = new TFile(argv[3]);
  TTree * Real_Tree = (TTree*)input_file_Real_DATA->Get("T");
// Get Information from MF and SRC
  TFile * input_file_mean = new TFile(argv[1], "update");   // input MF data
  TFile * input_file_SRC = new TFile(argv[2], "update");    // input SRC data
// Create output Trees/file for new weights
  TFile * output_weights = new TFile(argv[4], "RECREATE");                    // outfile of new weights
  TTree * outtree_weights = new TTree("updated_weights", "Weights");          // Will add MF and SRC trees as "friends"
  TTree * MF_Tree_Weights = new TTree("updated_MF_weights", "MF Weights");    // Holder to store MF weights
  TTree * SRC_Tree_Weights = new TTree("updated_SRC_weights", "SRC Weights"); // Holder to store SRC weights
    
// Get Trees for MF and SRC -> with approprate naming depending of the file
  TTree * Mean_Tree, * SRC_Tree;
  if (applied_simulator){
    Mean_Tree = (TTree*)input_file_mean->Get("T");
    SRC_Tree = (TTree*)input_file_SRC->Get("T");}
  else{
    Mean_Tree = (TTree*)input_file_mean->Get("genT");
    SRC_Tree = (TTree*)input_file_SRC->Get("genT");}

  
// Define Variables needed for Trees
  double Pbz = 4.461;      // Beam Energy GeV
  const int maxPart = 50.; // theoretical max length of momentum arrays (defined below)
  double pLead[3], pe[3], mom_x[maxPart], mom_y[maxPart], mom_z[maxPart];
  double updated_weight_SRC, updated_weight_MF, weight, Xb, Q2;
  int nParticles, Part_type[maxPart], lead_type;
  
// Get Data from Real Data Set
  Real_Tree->SetBranchAddress("mom_x"     , mom_x);        //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("mom_y"     , mom_y);        //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("mom_z"     , mom_z);        //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("nParticles", &nParticles);  //number of particles detected, integer
  Real_Tree->SetBranchAddress("Part_type" ,  Part_type );  //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
  Real_Tree->SetBranchAddress("Xb", &Xb);                  // Bjorken X
  Real_Tree->SetBranchAddress("Q2", &Q2);                  // Q_Squared

// Get Data from Mean Tree
  MF_Tree_Weights->Branch("updated_weight_MF",&updated_weight_MF,"updated_weight_MF/D");
  if (applied_simulator){
    Mean_Tree->SetBranchAddress("Part_type" , Part_type); //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
    Mean_Tree->SetBranchAddress("weighted", &weight);
    Mean_Tree->SetBranchAddress("mom_x"     , mom_x);     // momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("mom_y"     , mom_y);     // momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("mom_z"     , mom_z);     // momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("Q2", &Q2);               // Q_Squared
    Mean_Tree->SetBranchAddress("Xb", &Xb);}              // Bjorken Scaling Parameter
  else{
    Mean_Tree->SetBranchAddress("lead_type", &lead_type); // Code for Lead Nucleon
    Mean_Tree->SetBranchAddress("weight",&weight);
    Mean_Tree->SetBranchAddress("pLead", pLead);          // momentum vector of ejected nucleon (initial nucleon)
    Mean_Tree->SetBranchAddress("pe", pe);                // momentum vector of scattered electron (final)
    Mean_Tree->SetBranchAddress("QSq", &Q2);              // Q_Squared
    Mean_Tree->SetBranchAddress("xB", &Xb);}              // Bjorken Scaling Parameter
  
// Get Data from SRC Tree
  SRC_Tree_Weights->Branch("updated_weight_SRC",&updated_weight_SRC,"updated_weight_SRC/D");  // new weight to store
  if (applied_simulator){
    SRC_Tree->SetBranchAddress("Part_type" , Part_type);  //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
    SRC_Tree->SetBranchAddress("weighted", &weight);
    SRC_Tree->SetBranchAddress("mom_x"     , mom_x);      // momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("mom_y"     , mom_y);      // momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("mom_z"     , mom_z);      // momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("Q2", &Q2);                // Q_Squared
    SRC_Tree->SetBranchAddress("Xb", &Xb);}               // Bjorken Scaling Parameter
  else{
    SRC_Tree->SetBranchAddress("lead_type", &lead_type);  // Code for Lead Nucleon
    SRC_Tree->SetBranchAddress("weight",&weight);
    SRC_Tree->SetBranchAddress("pLead", pLead);           // momentum vector of ejected nucleon (1rst nucleon in SRC pair)
    SRC_Tree->SetBranchAddress("pe", pe);                 // momentum vector of scattered electron (final)
    SRC_Tree->SetBranchAddress("QSq", &Q2);               // Q_Squared
    SRC_Tree->SetBranchAddress("xB", &Xb);}               // Bjorken Scaling Parameter



  
// End Initialization of Files/Trees
  // ------------------------------------------------------------------------------------------------------------------- //
// Begin Program


// Look over small momentum segments to find optimal weighted scale at different points
// Scale should be in a well known region as it will be extended to whole data set
  std:: cout << "Looking at Small Momentum Segment Between: [" << start_p_scale << " , " << end_p_scale << "] MeV\n";
  
  
// Define Variables for Loop
  TVector3 P_miss_MF, P_miss_SRC;
  TTree *tree_array[2];
  tree_array[0] = Mean_Tree; tree_array[1] = SRC_Tree;
  double weight_sum[2], new_weight, scale;;
  weight_sum[0] = weight_sum[1] = 0.;
  double Scale;

// Loop over to get Real DATA Entries, apply cuts to make physical data -> NOTE: CUTS MUST MATCH IN FIND_MASS_COMBINED.root
// Making histograms of missing momentum to compare to SRC/MF data
  double Real_DATA_Entries = 0.;
  for(int i = 0; i < Real_Tree->GetEntries(); i++){
    Real_Tree->GetEntry(i);
    if (nParticles > 3. or nParticles < 2.) continue;
    if (Xb > xB_max_cut or Xb < xB_min_cut) continue;
    if (Q2 > QSq_max_cut or Q2 < QSq_min_cut) continue;
    if (Part_type[1] != 2212.) continue;
    TVector3 P_miss_real(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));
    P_miss_REAL_hist->Fill(P_miss_real.Mag()*1000.);
    // Only count the entry for scaling if it is in our region
    if (P_miss_real.Mag()*1000 < start_p_scale or P_miss_real.Mag()*1000 > end_p_scale) continue;
    Real_DATA_Entries +=1.;}

   
// Loop over all entries to find the total sum of the weight. Not saving the values yet
  for (int tree_type = 0; tree_type < 2; tree_type++){
    std::cout << "looping to find sum of weight \n";
    double total_entries = tree_array[tree_type]->GetEntries();
  for (int i = 0; i < total_entries; i++){
    tree_array[tree_type]->GetEvent(i);
    // Get missing momentum (MF label is arbitrary, for both SRC/MF)
    // Adding Same cuts as to Real DATA
    if (Xb > xB_max_cut or Xb < xB_min_cut) continue;
    if (Q2 > QSq_max_cut or Q2 < QSq_min_cut) continue;
    if (applied_simulator){
      if (Part_type[1] != 2212.) continue;
      P_miss_MF.SetXYZ(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));}
    else{
      if (lead_type != 2212.) continue;
      P_miss_MF.SetXYZ(pLead[0] + pe[0], pLead[1] + pe[1], pLead[2] - (Pbz - pe[2]));}
    if (P_miss_MF.Mag()*1000 < start_p_scale or P_miss_MF.Mag()*1000 > end_p_scale) continue;
    new_weight = weight/total_entries;
    weight_sum[tree_type] = weight_sum[tree_type] + new_weight;
  }
}

// Finding Scale Factor according to:
// total weight = Real_DATA_Entries = scale * (weight_MF + weight_SRC)
  Scale = Real_DATA_Entries/(weight_sum[0] + weight_sum[1]);
  std::cout << "Scale for Region: " << Scale << "\n";

// Loop over all entries in all trees to store final weight with scale (We will extend the scale to whole graph)
// Create graph of missing momentum distribution to compare to real data
  for (int tree_type = 0.; tree_type < 2.; tree_type++){
    std::cout << "looping to add to tree \n";
    double total_entries = tree_array[tree_type]->GetEntries();
  for (int i = 0; i < total_entries; i++){   
    tree_array[tree_type]->GetEvent(i);
    // Adding Same cuts as to Real DATA
    if (Xb > xB_max_cut or Xb < xB_min_cut) continue;
    if (Q2 > QSq_max_cut or Q2 < QSq_min_cut) continue;
    // Mean Field Route
    if (tree_type == 0.){
      // Get missing momentum
      if (applied_simulator){
	if (Part_type[1] != 2212.) continue;
	P_miss_MF.SetXYZ(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));}
      else{
	if (lead_type != 2212.) continue;
	P_miss_MF.SetXYZ(pLead[0] + pe[0], pLead[1] + pe[1], pLead[2] - (Pbz - pe[2]));}
      P_miss_MF_hist->Fill(P_miss_MF.Mag()*1000., weight*Scale/total_entries);
      updated_weight_MF = weight*Scale/total_entries;
      MF_Tree_Weights->Fill();}
    // SRC Route
    if (tree_type == 1.){
      // Get missing momentum
      if (applied_simulator){
	if (Part_type[1] != 2212.) continue;
	P_miss_SRC.SetXYZ(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));}
      else{
	if (lead_type != 2212.) continue;
	P_miss_SRC.SetXYZ(pLead[0] + pe[0], pLead[1] + pe[1], pLead[2] - (Pbz - pe[2]));}
      P_miss_SRC_hist->Fill(P_miss_SRC.Mag()*1000., weight*Scale/total_entries);
      updated_weight_SRC = weight*Scale/total_entries;
      SRC_Tree_Weights->Fill();}
  }
}
  

  
// Save it all in one Tree as a "friend"
  outtree_weights->AddFriend(SRC_Tree_Weights);
  outtree_weights->AddFriend(MF_Tree_Weights);
  
// Clean up
  // Write trees
  MF_Tree_Weights->Write();
  SRC_Tree_Weights->Write();
  outtree_weights->Write();
  // Write Plots
  P_miss_REAL_hist->Write();
  P_miss_SRC_hist->Write();
  P_miss_MF_hist->Write();
  // Close DATA
  input_file_Real_DATA->Close();
  input_file_mean->Close();
  input_file_SRC->Close();
  output_weights->Close();
  return 0.;
}
