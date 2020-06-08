#include <iostream>
#include <cmath>


#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>

#include "TString.h"

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"





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
  // Full P_miss Distribution
  const int total_bins = 100.;    // Total bins for P_miss histograms
  const int lower_bound = 0.;
  const int upper_bound = 1000.;
  // Ratio Graph Ranges
  const int section_width = 20;
  const int total_sections = (upper_bound-lower_bound)/section_width;
  // P_miss subsections to find weight
  const int start_MF = 0.;      // Initial P_miss value to consider in scaling (DATA comparison between real data)
  const int end_MF = 200.;      // Final P_miss calue to consider in scaling (DATA comparison between real data)
  const int start_SRC = 450.;   // Initial P_miss value to consider in scaling (DATA comparison between real data)
  const int end_SRC = 700.;     // Final P_miss calue to consider in scaling (DATA comparison between real data)
  // Cuts
  const bool cut_data = true;
  const double xB_min_cut = 1.15; // xB cut for data
  
// Make Histograms of Missing Momentum to Compare
  // Real
  TH1D * P_miss_REAL_hist = new TH1D("P_miss From REAL DATA","counts; Missing Momentum [MeV]", total_bins,  lower_bound, upper_bound);
  // SRC
  TH1D * P_miss_SRC_hist = new TH1D("P_miss From SRC DATA","counts; Missing Momentum [MeV]", total_bins,  lower_bound, upper_bound);
  TH1D * P_miss_SRC_FINAL = new TH1D("P_miss From SRC DATA FINAL Weight","counts; Missing Momentum [MeV]", total_bins, lower_bound, upper_bound);
  // MF
  TH1D * P_miss_MF_hist = new TH1D("P_miss From MF DATA","counts; Missing Momentum [MeV]", total_bins,  lower_bound, upper_bound);
  TH1D * P_miss_MF_FINAL = new TH1D("P_miss From MF DATA FINAL Weight","counts; Missing Momentum [MeV]", total_bins, lower_bound, upper_bound);
  // Combined Simulator
  TH1D * P_miss_combined_sim = new TH1D("P_miss From MF+SRC DATA FINAL","counts; Missing Momentum [MeV]", total_bins, lower_bound, upper_bound);
// Errors
  // P_miss
  P_miss_REAL_hist->Sumw2();
  // MF
  P_miss_MF_hist->Sumw2();
  P_miss_MF_FINAL->Sumw2();
  // SRC
  P_miss_SRC_hist->Sumw2();
  P_miss_SRC_FINAL->Sumw2();
  // Combined
  P_miss_combined_sim->Sumw2();
  

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
  double updated_weight_SRC, updated_weight_MF, weight, Xb;
  int nParticles, Part_type[maxPart], lead_type;
  
// Get Data from Real Data Set
  Real_Tree->SetBranchAddress("mom_x"     , mom_x);        //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("mom_y"     , mom_y);        //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("mom_z"     , mom_z);        //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Real_Tree->SetBranchAddress("nParticles", &nParticles);  //number of particles detected, integer
  Real_Tree->SetBranchAddress("Part_type" ,  Part_type );  //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
  Real_Tree->SetBranchAddress("Xb", &Xb);                  // Bjorken X

// Get Data from Mean Tree
  MF_Tree_Weights->Branch("updated_weight_MF",&updated_weight_MF,"updated_weight_MF/D");
  if (applied_simulator){
    Mean_Tree->SetBranchAddress("Part_type" , Part_type); //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
    Mean_Tree->SetBranchAddress("weighted", &weight);
    Mean_Tree->SetBranchAddress("mom_x"     , mom_x);     // momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("mom_y"     , mom_y);     // momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("mom_z"     , mom_z);     // momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Mean_Tree->SetBranchAddress("Xb", &Xb);}              // Bjorken Scaling Parameter
  else{
    Mean_Tree->SetBranchAddress("lead_type", &lead_type); // Code for Lead Nucleon
    Mean_Tree->SetBranchAddress("weight",&weight);
    Mean_Tree->SetBranchAddress("pLead", pLead);          // momentum vector of ejected nucleon (initial nucleon)
    Mean_Tree->SetBranchAddress("pe", pe);                // momentum vector of scattered electron (final)
    Mean_Tree->SetBranchAddress("xB", &Xb);}              // Bjorken Scaling Parameter
  
// Get Data from SRC Tree
  SRC_Tree_Weights->Branch("updated_weight_SRC",&updated_weight_SRC,"updated_weight_SRC/D");  // new weight to store
  if (applied_simulator){
    SRC_Tree->SetBranchAddress("Part_type" , Part_type);  //neutron (2112); proton (2212); pion (0: 111 ; -: -211 ; +: 211); electron (-11)
    SRC_Tree->SetBranchAddress("weighted", &weight);
    SRC_Tree->SetBranchAddress("mom_x"     , mom_x);      // momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("mom_y"     , mom_y);      // momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("mom_z"     , mom_z);      // momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    SRC_Tree->SetBranchAddress("Xb", &Xb);}               // Bjorken Scaling Parameter
  else{
    SRC_Tree->SetBranchAddress("lead_type", &lead_type);  // Code for Lead Nucleon
    SRC_Tree->SetBranchAddress("weight",&weight);
    SRC_Tree->SetBranchAddress("pLead", pLead);           // momentum vector of ejected nucleon (1rst nucleon in SRC pair)
    SRC_Tree->SetBranchAddress("pe", pe);                 // momentum vector of scattered electron (final)
    SRC_Tree->SetBranchAddress("xB", &Xb);}               // Bjorken Scaling Parameter



  
// End Initialization of Files/Trees
  // ------------------------------------------------------------------------------------------------------------------- //
// Begin Program

  
// Define Variables for Loop
  TVector3 P1_prime_TVec; double Pprime_over_Q, theta_P1_prime_q;
  int real_events_SRC = 0; int real_events_MF = 0;
// Loop over to get Real DATA Entries, apply cuts to make physical data -> NOTE: CUTS MUST MATCH IN FIND_MASS_COMBINED.root
  for(int i = 0; i < Real_Tree->GetEntries(); i++){
    Real_Tree->GetEntry(i);
    // Cuts to ONLY Real Data
    if (nParticles > 3. or nParticles < 2.) continue;
    // Cuts to Real and Sim DATA
    if (cut_data){
      if (Xb < xB_min_cut) continue;
      if (Part_type[1] != 2212.) continue;
      TVector3 q_TVec(-mom_x[0], -mom_y[0], Pbz - mom_z[0]); //turn to TVector, transfered momentum
      P1_prime_TVec.SetXYZ(mom_x[1], mom_y[1], mom_z[1]);          //turn to TVector, final nucleon momentum
      Pprime_over_Q = P1_prime_TVec.Mag() / q_TVec.Mag();
      if ( Pprime_over_Q > 0.96 or Pprime_over_Q < 0.62) continue; // selection criterion; taofeng    
      theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);              //give angle between vectors, Radians    
      if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;}         // selection criterion; taofeng
    
    // Add to Histograms
    TVector3 P_miss_real(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));
    P_miss_REAL_hist->Fill(P_miss_real.Mag()*1000.);
    // Only count the entry for scaling if it is in our region
    if (P_miss_real.Mag()*1000 > start_MF and P_miss_real.Mag()*1000 < end_MF) {real_events_MF += 1;}
    else if (P_miss_real.Mag()*1000 > start_SRC and P_miss_real.Mag()*1000 < end_SRC) {real_events_SRC += 1;}
  }

  
// Define Variables for Loop
  TVector3 P_miss_sim;
  TTree *tree_array[2];
  tree_array[0] = Mean_Tree; tree_array[1] = SRC_Tree;
  // Variables to count num MF and SRC in sections
  long double num_MF[total_sections+1], num_SRC[total_sections+1];  // Add 1 to capture edge effects
  for (int i=0 ; i<total_sections ; i++) {num_MF[i] = 0.; num_SRC[i] = 0.;}
  
  
// Loop over all entries in all trees
  double sum_SRC_weight_SRC = 0; double sum_SRC_weight_MF = 0;
  double sum_MF_weight_SRC = 0; double sum_MF_weight_MF = 0;
  double MF_events = 0; double SRC_events = 0;
  for (int tree_type = 0; tree_type < 2; tree_type++){
    std::cout << "looping through events \n";
    double total_entries = tree_array[tree_type]->GetEntries();
  for (int i = 0; i < total_entries; i++){
    tree_array[tree_type]->GetEvent(i);
    
    // Apply Cuts
    if (cut_data){
      if (Xb < xB_min_cut) continue;
      if (Part_type[1] != 2212.) continue;
      TVector3 q_TVec(-mom_x[0], -mom_y[0], Pbz - mom_z[0]); //turn to TVector, transfered momentum
      P1_prime_TVec.SetXYZ(mom_x[1], mom_y[1], mom_z[1]);          //turn to TVector, final nucleon momentum
      Pprime_over_Q = P1_prime_TVec.Mag() / q_TVec.Mag();
      if ( Pprime_over_Q > 0.96 or Pprime_over_Q < 0.62) continue; // selection criterion; taofeng    
      theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);              //give angle between vectors, Radians    
      if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;}        // selection criterion; taofeng

    // Calculate P_miss (add lead proton cut as that is true in E2A data.
    if (applied_simulator){
      if (Part_type[1] != 2212.) continue;
      P_miss_sim.SetXYZ(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));}
    else{
      if (lead_type != 2212.) continue;
      P_miss_sim.SetXYZ(pLead[0] + pe[0], pLead[1] + pe[1], pLead[2] - (Pbz - pe[2]));}
    
    // Only count the entry for scaling if it is in our region
    if (tree_type == 0){
      MF_events += 1;
      P_miss_MF_hist->Fill(P_miss_sim.Mag()*1000., weight);
      // For hist subsections
      if (P_miss_sim.Mag()*1000 > start_MF and P_miss_sim.Mag()*1000 < end_MF) {sum_MF_weight_MF += weight;}
      if (P_miss_sim.Mag()*1000 > start_SRC and P_miss_sim.Mag()*1000 < end_SRC) {sum_MF_weight_SRC += weight;}
      // Cut so we dont go over total_sections for the num_MF/SRC counting
      if (P_miss_sim.Mag()*1000 > upper_bound or P_miss_sim.Mag()*1000 < lower_bound) continue;
      num_MF[(int)(P_miss_sim.Mag()*1000/section_width)] += 1.;}
    else if (tree_type == 1){
      SRC_events += 1;
      P_miss_SRC_hist->Fill(P_miss_sim.Mag()*1000., weight);
      // For hist subsections
      if (P_miss_sim.Mag()*1000 > start_MF and P_miss_sim.Mag()*1000 < end_MF) {sum_SRC_weight_MF += weight;}
      if (P_miss_sim.Mag()*1000 > start_SRC and P_miss_sim.Mag()*1000 < end_SRC) {sum_SRC_weight_SRC += weight;}
      // Cut so we dont go over total_sections for the num_MF/SRC counting
      if (P_miss_sim.Mag()*1000 > upper_bound or P_miss_sim.Mag()*1000 < lower_bound) continue;
      num_SRC[(int)(P_miss_sim.Mag()*1000/section_width)] += 1.;}
  }
}

  // Find and thenuse the SRC Scale to Calculate MF Scale
  double SRC_scale = real_events_SRC/(sum_SRC_weight_SRC + sum_MF_weight_SRC);
  double MF_scale = (real_events_MF - (SRC_scale*sum_SRC_weight_MF))/(sum_MF_weight_MF);
  std::cout << "MF Scale: " << MF_scale << "\t" << "SRC Scale: " << SRC_scale  << "\n";
    
  
// Save new Weight as Friend
for (int tree_type = 0; tree_type < 2; tree_type++){
    std::cout << "looping through events to save new weights\n";
    double total_entries = tree_array[tree_type]->GetEntries();
  for (int i = 0; i < total_entries; i++){
    tree_array[tree_type]->GetEvent(i);
    // Apply Cuts, find P_miss
    if(cut_data){
      if (Xb < xB_min_cut) continue;
      if (Part_type[1] != 2212.) continue;
      TVector3 q_TVec(-mom_x[0], -mom_y[0], Pbz - mom_z[0]); //turn to TVector, transfered momentum
      P1_prime_TVec.SetXYZ(mom_x[1], mom_y[1], mom_z[1]);          //turn to TVector, final nucleon momentum
      Pprime_over_Q = P1_prime_TVec.Mag() / q_TVec.Mag();
      if ( Pprime_over_Q > 0.96 or Pprime_over_Q < 0.62) continue; // selection criterion; taofeng    
      theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);              //give angle between vectors, Radians    
      if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;}         // selection criterion; taofeng
    
    if (applied_simulator){
      if (Part_type[1] != 2212.) continue;
      P_miss_sim.SetXYZ(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));}
    else{
      if (lead_type != 2212.) continue;
      P_miss_sim.SetXYZ(pLead[0] + pe[0], pLead[1] + pe[1], pLead[2] - (Pbz - pe[2]));}
    
    // Only count the entry for scaling if it is in our region
    if (tree_type == 0){
      updated_weight_MF = weight*MF_scale;
      P_miss_combined_sim->Fill(P_miss_sim.Mag()*1000., updated_weight_MF);      
      P_miss_MF_FINAL->Fill(P_miss_sim.Mag()*1000., updated_weight_MF);      
      MF_Tree_Weights->Fill();}
    else if (tree_type == 1){
      updated_weight_SRC = weight*SRC_scale;
      P_miss_combined_sim->Fill(P_miss_sim.Mag()*1000., updated_weight_SRC);      
      P_miss_SRC_FINAL->Fill(P_miss_sim.Mag()*1000., updated_weight_SRC);
      SRC_Tree_Weights->Fill();}
   }
}


// Initialize variables
double P_miss_sect[total_sections];
double MF_SRC_ratio[total_sections];
double norm = SRC_events/MF_events;
 // Loop over all entries in all trees
  for (int tree_type = 0; tree_type < 2; tree_type++){
    std::cout << "looping through events \n";
  // Loop over and fill data
  for (int round = 0.; round < total_sections; round++){
    if ((int)num_MF[round] == 0.){
      MF_SRC_ratio[round] = 0.;}
    else if ((int)num_SRC[round] == 0.){
      std::cout << "No SRC Events at P_miss = : " << section_width*round + (section_width/2) + lower_bound << " -> Set ratio to 1000\n";
      MF_SRC_ratio[round] = 1000.;}
    else{
      MF_SRC_ratio[round] = norm*num_MF[round]/num_SRC[round];}
    P_miss_sect[round] = section_width*round + (section_width/2) + lower_bound;}
  }

// Graph Ratio of MN/SRC Number
  TCanvas *c2 = new TCanvas("Ratio of MF/SRC ", "Ratio Graph", 900, 900);
  TGraph *Ratio = new TGraph(total_sections, P_miss_sect, MF_SRC_ratio); // total points, x, y
  // Draw and write data
  Ratio->Draw("ALP");
  output_weights->cd();
  c2->Update();
  Ratio->Write();



// End of Progam
  // ------------------------------------------------------------------------------------------------------------------- //
// Output to user and clean up

  
// Save it all in one Tree as a "friend"
  outtree_weights->AddFriend(SRC_Tree_Weights);
  outtree_weights->AddFriend(MF_Tree_Weights);
  
// Clean up
  // Write trees
  MF_Tree_Weights->Write();
  SRC_Tree_Weights->Write();
  outtree_weights->Write();
  // Write histograms
  // P_miss
  P_miss_REAL_hist->Write();
  // MF
  P_miss_MF_hist->Write();
  P_miss_MF_FINAL->Write();
  // SRC
  P_miss_SRC_hist->Write();
  P_miss_SRC_FINAL->Write();
  // Combined
  P_miss_combined_sim->Write();
  // Close DATA
  input_file_Real_DATA->Close();
  input_file_mean->Close();
  input_file_SRC->Close();
  output_weights->Close();
// Close new weight file if added
   return 0;
}








/**
// Prefered method below, but not working.
R_Squared = corr.get_correlation(Real_HISTO, his_P1_before);
or
prob = hist1->KolmogorovTest(hist2,“N”);
**/






// Old Method
// -----------------------------------------------------------------------------------------//
// Find MF and SRC Weights
 
/**
// Desifne variables for Chi2 Test
  double MF_best_weight = -1.; double SRC_best_weight = -1;
  double Best_p_value = 0;
  double MF_weight, SRC_weight;
  const int max_rounds = 100;
  TH1F *proj_histo[(int) max_rounds];

  
// Create Histograms for every comparison
for (int round = 0.; round < max_rounds; round++){
  proj_histo[round] = NULL;
  std::stringstream histogramNameStream;
  histogramNameStream << "MF_weight_of_" << round/max_rounds << "_";
  std::string histogramName = histogramNameStream.str();
  proj_histo[round] = new TH1F( histogramName.c_str() , histogramName.c_str(),  total_bins, lower_bound, upper_bound);
 }

// Compare real data histograms to SRC + MF data
// Using different scaling weights in adding histograms to find best weight
//bool Best_Value_Found = false; // remains false if not good fit found


 TH1D *sim_histo[4];
 sim_histo[0] = P_miss_MF_hist_MF; sim_histo[1] = P_miss_SRC_hist_MF;
 sim_histo[2] = P_miss_MF_hist_SRC; sim_histo[3] = P_miss_SRC_hist_SRC;
 TH1D *Real_region[3];
 Real_region[0] = P_miss_REAL_hist_MF;
 Real_region[1] = P_miss_REAL_hist_SRC;

 
for (int tree_type = 3; tree_type > 0; tree_type -= 2){
  std::cout << "looping through events \n";
for (int round = 0.; round < max_rounds; round++){
  std::cout <<tree_type <<"\n";
  MF_weight = round/max_rounds;
  SRC_weight = 1 - MF_weight;
  
// Combine and Save Graphs with new weights:
   proj_histo[round]->Reset();
   proj_histo[round]->Add(sim_histo[tree_type], sim_histo[tree_type+1], MF_weight, SRC_weight); // MF data, SRC data, MF weight, SRC weight
   //proj_histo[round]->Write();  
   
// Compare histograms
  double res[total_bins];
  double p_value = Real_region[tree_type]->Chi2Test(proj_histo[round],"WW",res); // UW: unweight compared to weighted histogram (unweighted goes first). WW for both weighted. P for print.
  // If new fit is better, save the fit weight/statistics
  if (p_value > Best_p_value){
    //Best_Value_Found = true;
    Best_p_value = p_value;
    if (tree_type == 0){
      MF_best_weight = MF_weight;}
    else if (tree_type == 3){
      SRC_best_weight = SRC_weight;}}
 }
 }
**/

 
// End of Progam
  // ------------------------------------------------------------------------------------------------------------------- //
// Output to user and clean up
 
//std::cout << "Total SRC events: " << SRC_events << "\t Total Mean Field events: " << MF_events << "\n";
  //std::cout << "Best MF Scaling: " << MF_best_weight << "\t Best SRC Scaling: " << SRC_best_weight << "\t p-value of Fit: " << Best_p_value << "\n";
  //if (not Best_Value_Found){std:: cout << "Program did not work (or optimal MF scale is 0). No weighting was updated\n";}

