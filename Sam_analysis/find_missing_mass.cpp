#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"


// The purpose of this code is to calculate the total missing mass from mean_field and SRC pseudodata

int main(int argc, char** argv){
  bool new_weight, use_real_data;
  TFile * input_file_Real;
  TTree * Real_Tree;
  if (argc > 6. or argc < 5.){
    std::cout << "wrong number of arguments\n";
    std::cout << "\t Try ./find_missing_mass [input MF Data (root)] [input SRC Data (root)] [input Real Data (root) or type 1] [output file (root)] [Add Updated Weight Branch (Root file) or skip]\n";
    std::cout << "NOTE: MF and SRC input cannot be the same file. Change weights in code if this causes a problem\n";
    std::cout << "NOTE: If updated weight branch already added in previous code (as friend), you dont need to add it a second time here\n";
    return -1;
  }
  // Check if adding new weights
  if (argc == 5.){
    std::cout << "Using Weight in Input Root File. No Variables Changed\n";
    new_weight = false;}
  if (argc == 6.){
    std::cout << "Using New Weight in Updated Root File. No Other Variables Changed\n";
    new_weight = true;}
  // check if graphing real data
  if (atoi(argv[3]) == 1){
    std::cout << "Not graphing Real Data\n";
    use_real_data = false;}
  else{
    std::cout << "Graphing Real Data. If not ROOT file, error will be thrown later\n";
    use_real_data = true;}


// Create input and output files
  TFile * input_file_mean = new TFile(argv[1]);
  TFile * input_file_SRC = new TFile(argv[2]);
  TFile * output_file = new TFile(argv[4],"RECREATE");
  TFile * input_file_updated_weights;
  TTree * updated_weight_MF_tree, * updated_weight_SRC_tree;
  // add new weight file if inputted
  if (new_weight){
    input_file_updated_weights = new TFile(argv[5]);
    updated_weight_MF_tree = (TTree*)input_file_updated_weights->Get("updated_MF_weights");
    updated_weight_SRC_tree = (TTree*)input_file_updated_weights->Get("updated_SRC_weights");}
  // get real data file if added
  if (use_real_data){
    input_file_Real = new TFile(argv[3]);
    Real_Tree = (TTree*)input_file_Real->Get("T");}
  
//Make trees and histograms for the nuclei
  TTree * Mean_Tree = (TTree*)input_file_mean->Get("T");
  TTree * SRC_Tree = (TTree*)input_file_SRC->Get("T");


// Create Histograms of Missing Mass
  // Mass graph
  const int total_bins = 20.;    // Histogram bins. Later on, you want to make sure you can divide by section width
  const int M_miss_y_min = 0.;     
  const int M_miss_y_max = 2500.;
  const int M_miss_y_av_min = 860.;     
  const int M_miss_y_av_max = 1040.;
  const int P_miss_x_min = 0.;
  const int P_miss_x_max = 1000.; // make sure it can work with total sections
  // theta graph
  const int total_bins_theta = 20.;
  const int theta_min = 0.;     
  const int theta_max = 25;
  const int PQ_min = 0.;
  const int PQ_max = 1.; // make sure it can work with total sections
// Combined Contribution
  TH2D * his_P1_Mtf_av = new TH2D("P1_VS_Mtf_av","P1_VS_Mtf_mean;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_av_min, M_miss_y_av_max); // bin #, min1, min2, max
  TH2D * his_P1_Mtf = new TH2D("P1_VS_Mtf","P1_VS_Mtf;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_min, M_miss_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q = new TH2D("theta_VS_P1prime_q","theta_VS_P1prime_q;P1prime/q;theta", total_bins_theta, PQ_min, PQ_max, total_bins_theta, theta_min, theta_max); // bin #, min1, min2, max
// SRC contribution
  TH2D * his_P1_Mtf_SRC = new TH2D("P1_VS_Mtf_SRC","P1_VS_Mtf;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_min, M_miss_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q_SRC = new TH2D("theta_VS_P1prime_q_SRC","theta_VS_P1prime_q;P1prime/q;theta", total_bins_theta, PQ_min, PQ_max, total_bins_theta, theta_min, theta_max); // bin #, min1, min2, max
// Mean Contribution
  TH2D * his_P1_Mtf_MF = new TH2D("P1_VS_Mtf_MF","P1_VS_Mtf;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_min, M_miss_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q_MF = new TH2D("theta_VS_P1prime_q_MF","theta_VS_P1prime_q;P1prime/q;theta", total_bins_theta, PQ_min, PQ_max, total_bins_theta, theta_min, theta_max); // bin #, min1, min2, max
// Real Data Graphs for Comparison
  TH2D *his_P1_Mtf_REAL_av, *his_P1_Mtf_REAL;
  TH2D *his_theta_P1prime_q_REAL;
  if (use_real_data){
    his_P1_Mtf_REAL_av = new TH2D("P1_VS_Mtf_REAL_av","P1_VS_Mtf_mean;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_av_min, M_miss_y_av_max);
    his_P1_Mtf_REAL = new TH2D("P1_VS_Mtf_REAL","P1_VS_Mtf;P1;Mtf", total_bins, P_miss_x_min, P_miss_x_max, total_bins, M_miss_y_min, M_miss_y_max);
    his_theta_P1prime_q_REAL = new TH2D("theta_VS_P1prime_q_REAL","theta_VS_P1prime_q;P1prime/q;theta", total_bins_theta, PQ_min, PQ_max, total_bins_theta, theta_min, theta_max);}

// Sum of squares (error)
  // MTF
  his_P1_Mtf->Sumw2();
  his_P1_Mtf_av->Sumw2();
  his_P1_Mtf_SRC->Sumw2();
  his_P1_Mtf_MF->Sumw2();
  // Theta
  his_theta_P1prime_q->Sumw2();
  his_theta_P1prime_q_MF->Sumw2();
  his_theta_P1prime_q_SRC->Sumw2();
  // If using real
  if (use_real_data){
    his_P1_Mtf_REAL->Sumw2();
    his_P1_Mtf_REAL_av->Sumw2();
    his_theta_P1prime_q_REAL->Sumw2();}

  
// Define Variables from massT
  const double me = 0.000511;
  const double Pbz = 4.461; // THIS SHOULD BE BEAM ENERGY, GeV! CHANGE GENERATOR IF NOT
  const int maxPart = 50.;
  double mom_x[maxPart], mom_y[maxPart], mom_z[maxPart];
  int Part_type[maxPart];
  double Xb, Q2, weighted; 

//Variables we are taking from pseudo_skim_tree for Mean Field
  Mean_Tree->SetBranchAddress("Part_type" ,  Part_type ); //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  Mean_Tree->SetBranchAddress("mom_x"     ,  mom_x     ); //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("mom_y"     ,  mom_y     ); //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("mom_z"     ,  mom_z     ); //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("Q2", &Q2);  // Q_Squared
  Mean_Tree->SetBranchAddress("Xb"        , &Xb        ); // Bjorken X
//Variables we are taking for SRC
  SRC_Tree->SetBranchAddress("Part_type" ,  Part_type ); //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  SRC_Tree->SetBranchAddress("mom_x"     ,  mom_x     ); //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("mom_y"     ,  mom_y     ); //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("mom_z"     ,  mom_z     ); //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("Q2", &Q2);  // Q_Squared
  SRC_Tree->SetBranchAddress("Xb"        , &Xb        ); // Bjorken X
//Variables we are taking from Real tree is provided
  if (use_real_data){
    Real_Tree->SetBranchAddress("Part_type" ,  Part_type );  //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
    Real_Tree->SetBranchAddress("mom_x"     ,  mom_x     );  //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Real_Tree->SetBranchAddress("mom_y"     ,  mom_y     );  //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Real_Tree->SetBranchAddress("mom_z"     ,  mom_z     );  //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
    Real_Tree->SetBranchAddress("Q2", &Q2);  // Q_Squared
    Real_Tree->SetBranchAddress("Xb"        , &Xb        );} // Bjorken X
// If Using updated weights
  if (new_weight){
    // MF Weight
    updated_weight_MF_tree->SetBranchAddress("updated_weight_MF", &weighted);
    Mean_Tree->AddFriend(updated_weight_MF_tree);
    // SRC Weight
    updated_weight_SRC_tree->SetBranchAddress("updated_weight_SRC", &weighted);
    SRC_Tree->AddFriend(updated_weight_SRC_tree);}
  else{
    Mean_Tree->SetBranchAddress("weighted", &weighted);
    SRC_Tree->SetBranchAddress("weighted", &weighted);}

  
// End Initialization of Variables/Files/Trees
  // ------------------------------------------------------------------------------------------------------------------- //
// Begin Program

  
// Define Variables for Loop
  double Mn, theta_P1_prime_q, Pprime_over_Q;
  TVector3 P1_prime_TVec;
  // Create array of trees to loop over
  int tree_size = 2;
  TTree *tree_array[3];
  tree_array[0] = Mean_Tree;
  tree_array[1] = SRC_Tree;
  if (use_real_data){
    tree_size += 1.;
    tree_array[2] = Real_Tree;}

// Loop over all entries in all trees
  for (int tree_type = 0.; tree_type < tree_size; tree_type++){
    std::cout << "looping through events \n";
  for (int i = 0.; i < tree_array[tree_type]->GetEntries(); i++){   
    tree_array[tree_type]->GetEvent(i);

// Apply SRC-like Kinematic Selection Criteria from Taofeng     
    if (Xb < 1.15) continue;
    if (Q2 > 4.1 or Q2 < 0.5) continue;
// Find out which of the pair is viable; if both: take first viable proton (unsure how to deal with both now)
    bool nucleon_test = false;
    int recorded_nucleon_index;
    TVector3 q_TVec(-mom_x[0], -mom_y[0], Pbz - mom_z[0]); //turn to TVector, transfered momentum
    for (int i = 1; i < 2./**nParticles**/; i++){ // look at first ejected particle for now
      P1_prime_TVec.SetXYZ(mom_x[i], mom_y[i], mom_z[i]);          //turn to TVector, final nucleon momentum
      Pprime_over_Q = P1_prime_TVec.Mag() / q_TVec.Mag();
      if ( Pprime_over_Q > 0.96 or Pprime_over_Q < 0.62) continue; // selection criterion; taofeng    
      theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);              //give angle between vectors, Radians    
      if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;         // selection criterion; taofeng
      // If we get to this point than the nucleon matched our criterion, so we break out of the 4-loop
      // Won't test other nucleon for now
      nucleon_test = true;
      recorded_nucleon_index = i;
      break;}
    if (not nucleon_test) continue; // go onto next event if no nucleon matched criteria

    
// Label the mass of exiting nucleon and nucleus
    if (Part_type[recorded_nucleon_index] == 2212.){
      Mn = 0.93827231;} // Mass of exiting proton; Gev
    else if (Part_type[recorded_nucleon_index] == 2112.){
      Mn = 0.93957;}    // Mass of exiting neutron; Gev
    else if (Part_type[recorded_nucleon_index] == -211.){
      Mn = 0.13957061;}    // Mass of exiting pion(-1); Gev
    else if (Part_type[recorded_nucleon_index] == 111){
      Mn = 0.134977;}    // Mass of exiting pion(0); Gev
    else if (Part_type[recorded_nucleon_index] == 211.){
      Mn = 0.13957061;}    // Mass of exiting pion(1); Gev **/
    else{
      std::cout << "no mass; Nucleon Type Not in Database \n";
      continue;}

    
// Define Variables in Lorentz Vector to calculate missing mass
    double Eb = sqrt(Pbz*Pbz + me*me);  // initial electron beam energy (about = P_z); Gev
    TVector3 Eb_prime_TVec(mom_x[0], mom_y[0], mom_z[0]); //turn to TVector, transfered momentum
    double Eb_prime = sqrt(Eb_prime_TVec.Mag2() + me*me);
    double w = Eb - Eb_prime;
    TLorentzVector P1_prime_final;
    P1_prime_final.SetPxPyPzE(P1_prime_TVec[0], P1_prime_TVec[1], P1_prime_TVec[2], sqrt(P1_prime_TVec.Mag2() + Mn*Mn)); //Px, Py, Pz, E
    TLorentzVector q_photon;
    q_photon.SetPxPyPzE(q_TVec[0], q_TVec[1], q_TVec[2], w); //Px, Py, Pz, E
    TLorentzVector P1_initial_pair;
    P1_initial_pair.SetPxPyPzE(0, 0, 0, 2*Mn); //Px, Py, Pz, E
    // Calculate Missing Mass
    double Missing_Mass = (q_photon + P1_initial_pair - P1_prime_final).Mag();


// Fill in Historams with Data    
    TVector3 P_miss = P1_prime_TVec - q_TVec; // Should = mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0])
    if ( tree_type == 0. ){
      his_P1_Mtf->Fill(P_miss.Mag()*1000., Missing_Mass*1000., weighted);
      his_P1_Mtf_MF->Fill(P_miss.Mag()*1000., Missing_Mass*1000., weighted);
      his_theta_P1prime_q->Fill(Pprime_over_Q, theta_P1_prime_q*(180/M_PI), weighted);
      his_theta_P1prime_q_MF->Fill(Pprime_over_Q, theta_P1_prime_q*(180/M_PI), weighted);}
    else if ( tree_type == 1. ){
      his_P1_Mtf->Fill(P_miss.Mag()*1000., Missing_Mass*1000., weighted);
      his_P1_Mtf_SRC->Fill(P_miss.Mag()*1000., Missing_Mass*1000., weighted);
      his_theta_P1prime_q->Fill(Pprime_over_Q, theta_P1_prime_q*(180/M_PI), weighted);
      his_theta_P1prime_q_SRC->Fill(Pprime_over_Q, theta_P1_prime_q*(180/M_PI), weighted);}
    else if ( tree_type == 2. ){
      weighted = 1.; // No weighted for Real Data;
      his_P1_Mtf_REAL->Fill(P_miss.Mag()*1000., Missing_Mass*1000., weighted);
      his_theta_P1prime_q_REAL->Fill(Pprime_over_Q, theta_P1_prime_q*(180/M_PI), weighted);}
 }}


  
// End of Data Mining
  // ------------------------------------------------------------------------------------------------------------------- //
// Calcuate the Average of Missing Mass Data for some P_miss region
  
  
// Find mean for every x point and replot
  int section_width = 50.;
  int total_sections = P_miss_x_max/section_width;
  double bin_per_section = total_bins/total_sections; // check to make sure divides
  double x[(int) total_sections];
  double y[(int) total_sections];
  double ex[(int) total_sections];
  double ey[(int) total_sections];
  TH1F *proj_histo[total_sections];
  for (int i=0 ; i<total_sections ; i++) proj_histo[i] = NULL;

// Create Histograms for every section (projection section) of the graph
for (int round = 0.; round < total_sections; round++){
  std::stringstream histogramNameStream;
  histogramNameStream << "Projection_from_" << (section_width*round) << "_to_" << (section_width*(round+1.));
  std::string histogramName = histogramNameStream.str();
  proj_histo[round] = new TH1F( histogramName.c_str() , ":(" ,  total_bins, M_miss_y_min, M_miss_y_max);
 }


// Variables for loop
 TH2D *M_miss_hist_array[2];
 M_miss_hist_array[0] = his_P1_Mtf;
 TH2D *M_miss_av_hist_array[2];
 M_miss_av_hist_array[0] = his_P1_Mtf_av;
 TGraphAsymmErrors *Mtf_error_array[2];
 TCanvas *canvas_array1[2], *canvas_array2[2];
 TString Graph_Names[2];
 Graph_Names[0] = "Mean Field + SRC Missing Mass";
 Graph_Names[1] = "Real Data Missing Mass";
 int array_length = 1.;
 // Add in loop if real data present
 if (use_real_data){
   array_length = 2.;
   M_miss_hist_array[1] = his_P1_Mtf_REAL;
   M_miss_av_hist_array[1] = his_P1_Mtf_REAL_av;}
// Loop over all entries to find the total sum of the weight. Not saving the values yet
for (int hist_type = 0.; hist_type < array_length; hist_type++){
  std::cout << "looping to find projected missing mass average \n";
    
// Recreate graph with mean plotted
std::ostringstream Canvas_Name1;
Canvas_Name1 << hist_type;
canvas_array1[hist_type] = new TCanvas(Canvas_Name1.str().c_str(),"c1",900,900);// Mass_x_min, Mass_x_max, Mass_y_min, Mass_y_max);
for (int round = 0.; round < total_sections; round++){
  // project the Y axis
  TH1D * proj_Mtf = M_miss_hist_array[hist_type]->ProjectionY(":)", bin_per_section*round + 1, bin_per_section*(round+1.)); // name, first bin, last bin

  // fit the function and save for reference
  //proj_histo[round] = NULL;
  proj_histo[round]->Add(proj_Mtf);
  proj_histo[round]->Fit("gaus", "QEM"); // fit the function to a gaussian
  // Write each fit
  if (hist_type == 0. and not new_weight and not use_real_data){
    proj_histo[round]->Write();}
  
  // Get values for Mass Graph
  double Mtf_bin_mean = proj_histo[round]->GetMean(); // find mean given in bin graph
  double Mtf_bin_std = proj_histo[round]->GetStdDev();
  int num_points = proj_histo[round]->GetEntries(); // gives number of entries as double -> turn to int
  x[round] = (round+0.5)*section_width;  
  y[round] = Mtf_bin_mean;
  ex[round] = 0.;
  ey[round] = Mtf_bin_std/sqrt(num_points - 1);
  M_miss_av_hist_array[hist_type]->Fill(x[round], y[round]);

  // Pring out fit for combined data
  if (hist_type == 0){
    std::cout << "MF/SRC: P1 Range: [" << section_width*round << ", " << section_width*(round+1.) << "]          "
	      << "fit std:  " <<  Mtf_bin_std << "\t fit mean:  " <<  Mtf_bin_mean
	      << "\t number of points in sample:  " << num_points << "\n";}
  proj_Mtf->Reset("ICESM");
  proj_histo[round]->Reset("ICESM");
 }
 
// Make mean graph pretty
  M_miss_av_hist_array[hist_type]->SetMarkerStyle(8);  
  M_miss_av_hist_array[hist_type]->SetMarkerSize(1);   
  M_miss_av_hist_array[hist_type]->SetMarkerColor(4);  
 
// Error Bar Graph of Mtf Mean
  std::ostringstream Canvas_Name;
  Canvas_Name << Graph_Names[hist_type];
  canvas_array2[hist_type] = new TCanvas(Canvas_Name.str().c_str(), Canvas_Name.str().c_str(), 900, 900);
  Mtf_error_array[hist_type] = new TGraphAsymmErrors(total_sections, x, y, ex, ex, ey, ey);
// Make Pretty and Draw   
   Mtf_error_array[hist_type]->SetTitle(Canvas_Name.str().c_str());
   Mtf_error_array[hist_type]->SetName(Canvas_Name.str().c_str());
   Mtf_error_array[hist_type]->SetLineColor(4);
   Mtf_error_array[hist_type]->SetMarkerStyle(8);
   Mtf_error_array[hist_type]->SetMarkerSize(1);                               
   Mtf_error_array[hist_type]->SetMarkerColor(4);
   Mtf_error_array[hist_type]->GetYaxis()->SetTitle("Missing Mass [MeV]");
   Mtf_error_array[hist_type]->GetYaxis()->SetLabelSize(0.03);     
   Mtf_error_array[hist_type]->GetYaxis()->SetTitleOffset(1.3);                
   Mtf_error_array[hist_type]->GetYaxis()->CenterTitle(true);                  
   Mtf_error_array[hist_type]->GetXaxis()->SetTitle("Missing Momentum [MeV]"); 
   Mtf_error_array[hist_type]->GetXaxis()->SetLabelSize(0.03);     
   Mtf_error_array[hist_type]->GetXaxis()->CenterTitle(true);                  
   Mtf_error_array[hist_type]->SetMinimum(860);                                
   Mtf_error_array[hist_type]->SetMaximum(1040);                               
   Mtf_error_array[hist_type]->SetFillColor(6);
   Mtf_error_array[hist_type]->SetFillStyle(3005);
   Mtf_error_array[hist_type]->Draw("ALP");
   // Write Missing Mass Average Plots with Error Bars
   output_file->cd();
   canvas_array2[hist_type]->Update();
   Mtf_error_array[hist_type]->Write();
   
}

// End of Missing Mass Projects
  // ------------------------------------------------------------------------------------------------------------------- //
// Clean Up Files
  
// Clean Up Files
    // Write Missing Mass Scatter Plots
    his_P1_Mtf->Write();
    his_P1_Mtf_MF->Write();
    his_P1_Mtf_SRC->Write();
    // Write Missing Mass Average Plots
    his_P1_Mtf_av->Write();
    // Write theta P'/q VS q plots  
    his_theta_P1prime_q_SRC->Write();
    his_theta_P1prime_q_MF->Write();
    his_theta_P1prime_q->Write();
    if (use_real_data){
      his_theta_P1prime_q_REAL->Write();
      his_P1_Mtf_REAL_av->Write();
      his_P1_Mtf_REAL->Write();}
    
    output_file->Close();
    // Close Input Files
    input_file_mean->Close();
    input_file_SRC->Close();
    if (use_real_data){
      input_file_Real->Close();}
    if (new_weight){
      input_file_updated_weights->Close();}
    return 0;
}

// End Program
  // ------------------------------------------------------------------------------------------------------------------- //
// Optional Additions

/**
// Add Gaussian Fit:
  #include "TF1.h"
  double fit_Mtf_bin_mean = 0;
  double fit_Mtf_bin_std = 0.;
  TF1 *fit = (TF1 *)proj_histo[round]->GetFunction("gaus");
  if(fit != NULL){
    fit_Mtf_bin_mean = fit->GetParameter(1);
    fit_Mtf_bin_std = fit->GetParameter(2);}
  else{
    std:: cout << "No Data to Fit. Set: Mean = STD = 0 \n";} **/
