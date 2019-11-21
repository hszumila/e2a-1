#include <iostream>
#include <cmath>

#include <string>
#include <vector>
#include <math.h>

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
  if (argc != 4){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try ./find_mass_combined [input MF Data = DATA/simulator_MF.root] [input SRC Data = DATA/simulator_SRC.root] [output = DATA/find_mass_combined.root]\n";
    std::cout << "NOTE: MF and SRC input cannot be the same file. Change weights in code if this causes a problem\n";
    return -1;
  }

  
// Create input and output files
  TFile * input_file_mean = new TFile(argv[1]);
  TFile * input_file_SRC = new TFile(argv[2]);
  TFile * output_file = new TFile(argv[3],"RECREATE");

// Label constant variables to use in program:
  const double total_bins = 20.; // Histogram bins. Later on, you want to make sure you can divide by section width
  const int Mass_y_min = 860.;     
  const int Mass_y_max = 1040.;
  const int Mass_x_min = 0.;
  const int Mass_x_max = 1000.; // make sure it can work with total sections
  const double mean_field_weight = 0.8;
  const double SRC_weight = 0.2;
  const double tree_weight[2] = { mean_field_weight, SRC_weight };
  
//Make trees and histograms for the nuclei
  TTree * Mean_Tree = (TTree*)input_file_mean->Get("T");
  TTree * SRC_Tree = (TTree*)input_file_SRC->Get("T");
  // Combined Contribution
  TH2D * his_P1_Mtf_av = new TH2D("P1_VS_Mtf_av","P1_VS_Mtf_mean;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_P1_Mtf = new TH2D("P1_VS_Mtf","P1_VS_Mtf;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q = new TH2D("theta_VS_P1prime_q","theta_VS_P1prime_q;P1prime/q;theta", 1000, 0, 1.5, 1000, 0, 200); // bin #, min1, min2, max
  TH1D * his_Q2_weight = new TH1D("Q2_weight","Q2 [GeV^2];Counts", 20, 0., 4.1);
  TH1D * his_Xb = new TH1D("Xb","Xb;Counts", 20, 1., 2.);
  TH1D * his_P1_before = new TH1D("P1 before: no cuts","P1;Counts", 20, 0., 1000.);
  TH1D * his_P1_after = new TH1D("P1 after: with cuts","P1;Counts", 20, 0., 1000.);
  // SRC contribution
  TH2D * his_P1_Mtf_SRC = new TH2D("P1_VS_Mtf_SRC","P1_VS_Mtf;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q_SRC = new TH2D("theta_VS_P1prime_q_SRC","theta_VS_P1prime_q;P1prime/q;theta", 1000, 0, 1.5, 1000, 0, 200); // bin #, min1, min2, max
  TH1D * his_Q2_weight_SRC = new TH1D("Q2_weight_SRC","Q2 [GeV^2];Counts", 20, 0., 4.1);
  TH1D * his_Xb_SRC = new TH1D("Xb: SRC","Xb;Counts", 20, 1., 2.);
  TH1D * his_P1_before_SRC = new TH1D("P1 before, SRC: no cuts","P1;Counts", 20, 0., 1000.);
  TH1D * his_P1_after_SRC = new TH1D("P1 after, SRC: with cuts","P1;Counts", 20, 0., 1000.);
  // Mean Contribution
  TH2D * his_P1_Mtf_MF = new TH2D("P1_VS_Mtf_mean","P1_VS_Mtf;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q_MF = new TH2D("theta_VS_P1prime_q_mean","theta_VS_P1prime_q;P1prime/q;theta", 1000, 0, 1.5, 1000, 0, 200); // bin #, min1, min2, max
  TH1D * his_Q2_weight_MF = new TH1D("Q2_weight_mean","Q2 [GeV^2];Counts", 20, 0., 4.1);
  TH1D * his_Xb_MF = new TH1D("Xb: MF","Xb;Counts", 20, 1., 2.);
  TH1D * his_P1_before_MF = new TH1D("P1 before, MF: no cuts","P1;Counts", 20, 0., 1000.);
  TH1D * his_P1_after_MF = new TH1D("P1 after, MF: with cuts","P1;Counts", 20, 0., 1000.);

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
  // Q2
  his_Q2_weight->Sumw2();  
  his_Q2_weight_SRC->Sumw2(); 
  his_Q2_weight_MF->Sumw2();
  // Xb
  his_Xb->Sumw2();
  his_Xb_SRC->Sumw2();
  his_Xb_MF ->Sumw2();
  // P1
  his_P1_before->Sumw2();
  his_P1_after->Sumw2();
  his_P1_before_SRC->Sumw2();
  his_P1_after_SRC->Sumw2();
  his_P1_before_MF->Sumw2();
  his_P1_after_MF->Sumw2();

  
// Define Variables from massT
  double Pbz = 4.461; // THIS SHOULD BE BEAM ENERGY, GeV! CHANGE GENERATOR IF NOT
  const int maxPart = 50.;
  double mom_x[maxPart], mom_y[maxPart], mom_z[maxPart], vtx_z_cor[maxPart];
  int nParticles, Part_type[maxPart];
  double Xb, Q2; 

//Variables we are taking from pseudo_skim_tree for Mean Field
  Mean_Tree->SetBranchAddress("nParticles", &nParticles); //number of particles detected, integer - for me it will always be two  
  Mean_Tree->SetBranchAddress("Part_type" ,  Part_type ); //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  Mean_Tree->SetBranchAddress("vtx_z_cor" ,  vtx_z_cor ); //vertext Z corrected, arrays of integer, length of nParticles
  Mean_Tree->SetBranchAddress("mom_x"     ,  mom_x     ); //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("mom_y"     ,  mom_y     ); //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("mom_z"     ,  mom_z     ); //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  Mean_Tree->SetBranchAddress("Xb"        , &Xb        ); // Bjorken X
  Mean_Tree->SetBranchAddress("Q2"        , &Q2        ); // Q_2
//Variables we are taking for SRC
  SRC_Tree->SetBranchAddress("nParticles", &nParticles); //number of particles detected, integer - for me it will always be two  
  SRC_Tree->SetBranchAddress("Part_type" ,  Part_type ); //neutron (2112) or proton (2212) or pion (0: 111 ; -: -211 ; +: 211) or electron (-11)
  SRC_Tree->SetBranchAddress("vtx_z_cor" ,  vtx_z_cor ); //vertext Z corrected, arrays of integer, length of nParticles
  SRC_Tree->SetBranchAddress("mom_x"     ,  mom_x     ); //momentum in x direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("mom_y"     ,  mom_y     ); //momentum in z direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("mom_z"     ,  mom_z     ); //momentum in y direction, arrays of double, length of nParticles, GeV, final nucleon momentum
  SRC_Tree->SetBranchAddress("Xb"        , &Xb        ); // Bjorken X
  SRC_Tree->SetBranchAddress("Q2"        , &Q2        ); // Q_2

  
// Define Variables for Loop
  double Mn, theta_P1_prime_q, weighted;
  TVector3 P1_prime_TVec;
  TTree *tree_array[2];
  tree_array[0] = Mean_Tree;
  tree_array[1] = SRC_Tree;
  double SRC_events, MF_events = 0;

// Loop over all entries in all trees
  for (int tree_type = 0; tree_type < 2; tree_type++){
    std::cout << "looping through events \n";
  for (int i = 0; i < tree_array[tree_type]->GetEntries(); i++){   
    tree_array[tree_type]->GetEvent(i);

    // delete later
    TVector3 P1_holder(mom_x[1] + mom_x[0], mom_y[1] + mom_y[0], mom_z[1] - (Pbz - mom_z[0]));
    if ( tree_type == 0. ){
      his_P1_before_MF->Fill(P1_holder.Mag()*1000, weighted * tree_weight[tree_type]);}
    else{
      his_P1_before_SRC->Fill(P1_holder.Mag()*1000, weighted * tree_weight[tree_type]);}
    // end delete later
      
    if (Xb < 1.15) continue;
    TVector3 q_TVec(-mom_x[0], -mom_y[0], Pbz - mom_z[0]); //turn to TVector, transfered momentum
    double w = sqrt(q_TVec.Mag2() - Q2);
    
// Find out which of the pair is viable; if both: take first viable one (unsure how to deal with both now)
    bool nucleon_test = false;
    int recorded_nucleon;
    for (int i = 1; i < 2./**nParticles**/; i++){ // look at ejected particles for one event for now
      P1_prime_TVec.SetX(mom_x[i]); P1_prime_TVec.SetY(mom_y[i]); P1_prime_TVec.SetZ(mom_z[i]);             //turn to TVector, final nucleon momentum
      if ((P1_prime_TVec.Mag()/q_TVec.Mag()) > 0.96 or (P1_prime_TVec.Mag()/q_TVec.Mag()) < 0.62) continue; // selection criterion; taofeng    
      theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);                                                       //give angle between vectors, Radians    
      if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;                                                  // selection criterion; taofeng
      // If we get to this point than the nucleon matched our criterion, so we break out of the 4-loop
      // Won't test other nucleon for now
      nucleon_test = true;
      recorded_nucleon = i;
      break;}
    if (not nucleon_test) continue; // go onto next event if no nucleon matched criteria

// Label the mass of exiting nucleon and nucleus
    if (Part_type[recorded_nucleon] == 2212.){
      Mn = 0.93827231;} // Mass of exiting proton; Gev
    else if (Part_type[recorded_nucleon] == 2112.){
      Mn = 0.93957;}    // Mass of exiting neutron; Gev
    else if (Part_type[recorded_nucleon] == -211.){
      Mn = 0.13957061;}    // Mass of exiting pion(-1); Gev
    else if (Part_type[recorded_nucleon] == 111){
      Mn = 0.134977;}    // Mass of exiting pion(0); Gev
    else if (Part_type[recorded_nucleon] == 211.){
      Mn = 0.13957061;}    // Mass of exiting pion(1); Gev **/
    else{
      std::cout << "no mass; Nucleon Type Not in Database \n";
      continue;}

   
// Define Variables in Lorentz Vector
    TLorentzVector P1_prime_final;
    P1_prime_final.SetPxPyPzE(P1_prime_TVec[0], P1_prime_TVec[1], P1_prime_TVec[2], sqrt(P1_prime_TVec.Mag2() + Mn*Mn)); //Px, Py, Pz, E
    TLorentzVector q_photon;
    q_photon.SetPxPyPzE(q_TVec[0], q_TVec[1], q_TVec[2], w); //Px, Py, Pz, E
    TLorentzVector P1_initial_pair;
    P1_initial_pair.SetPxPyPzE(0, 0, 0, 2*Mn); //Px, Py, Pz, E
    // Calculate Missing Mass
    double Mtf = sqrt((q_photon + P1_initial_pair - P1_prime_final).Mag2());

    
// FInd weight for event
// Real Data doesnt have weight (weight = 1). Only get weight if none provided (simulated data)
  double weight, weighted;
  TBranch * weight_branch = tree_array[tree_type]->GetBranch("weighted");
  if (weight_branch){
    tree_array[tree_type]->SetBranchAddress("weighted"  , &weighted  );
    weight = weighted;}//
  else{
    std::cout << "no weight; Real DATA -> weight = 1 \n";
    weight = 1;}
  if (weight * tree_weight[tree_type] == 0.) continue;

// Create histograms
    TVector3 P1_TVec = P1_prime_TVec - q_TVec;
    if ( tree_type == 0. ){
      MF_events += 1.;
      his_P1_Mtf_MF->Fill(P1_TVec.Mag()*1000., Mtf*1000., weight * tree_weight[tree_type]);
      his_theta_P1prime_q_MF->Fill((P1_prime_TVec.Mag()/q_TVec.Mag()), theta_P1_prime_q*(180/M_PI), weight * tree_weight[tree_type]);
      his_Q2_weight_MF->Fill(Q2, weight * tree_weight[tree_type]);
      his_Xb_MF->Fill(Xb, weight * tree_weight[tree_type]);
      his_P1_after_MF->Fill(P1_TVec.Mag()*1000, weight * tree_weight[tree_type]);}
    else if ( tree_type == 1. ){
      SRC_events += 1.;
      his_P1_Mtf_SRC->Fill(P1_TVec.Mag()*1000., Mtf*1000., weight * tree_weight[tree_type]);
      his_theta_P1prime_q_SRC->Fill((P1_prime_TVec.Mag()/q_TVec.Mag()), theta_P1_prime_q*(180/M_PI), weight * tree_weight[tree_type]);
      his_Q2_weight_SRC->Fill(Q2, weight * tree_weight[tree_type]);
      his_Xb_SRC->Fill(Xb, weight * tree_weight[tree_type]);
      his_P1_after_SRC->Fill(P1_TVec.Mag()*1000, weight * tree_weight[tree_type]);}
  }}

// Combine Graphs with weights:
  // Add SRC (scaled up to MF number)
  his_P1_Mtf->Add(his_P1_Mtf_SRC, MF_events); // add histogram, weight those values
  his_theta_P1prime_q->Add(his_theta_P1prime_q_SRC, MF_events);
  his_Q2_weight->Add(his_Q2_weight_SRC, MF_events);
  his_Xb->Add(his_Xb_SRC, MF_events);
  his_P1_before->Add(his_P1_before_SRC, MF_events);
  his_P1_after->Add(his_P1_after_SRC, MF_events);

  // Add MF (scaled up to SRC number)
  his_P1_Mtf->Add(his_P1_Mtf_MF, SRC_events);
  his_theta_P1prime_q->Add(his_theta_P1prime_q_MF, SRC_events);
  his_Q2_weight->Add(his_Q2_weight_MF, SRC_events);
  his_Xb->Add(his_Xb_MF, SRC_events);
  his_P1_before->Add(his_P1_before_MF, SRC_events);
  his_P1_after->Add(his_P1_after_MF, SRC_events);
  
  std::cout << "SRC events: " << SRC_events << "     Mean Field events: " << MF_events << "\n";

  // ------------------------------------------------------------------------------------------------------------------- //

  
// Find mean for every x point and replot
  int section_width = 50.;
  int total_sections = Mass_x_max/section_width;
  double bin_per_section = total_bins/total_sections; // check to make sure divides
  double x[(int) total_sections];
  double y[(int) total_sections];
  double ex[(int) total_sections];
  double ey[(int) total_sections];
  TH1F *proj_histo[total_sections];
  for (int i=0 ; i<total_sections ; i++) proj_histo[i]=NULL;

// Create Histograms for every section (projection section) of the graph
for (int round = 0.; round < total_sections; round++){
  std::string histogramName = "projection_from_" + std::to_string(section_width*round) + "_to_" + std::to_string(section_width*(round+1.));
  proj_histo[round] = new TH1F( histogramName.c_str() , ":(" ,  total_bins, Mass_y_min, Mass_y_max);
 }

// Recreate graph with mean plotted
TCanvas *c1 = new TCanvas("c1","c1",900,900);// Mass_x_min, Mass_x_max, Mass_y_min, Mass_y_max);
for (int round = 0.; round < total_sections; round++){
  // project the Y axis
  TH1D * proj_Mtf = his_P1_Mtf->ProjectionY(":)", bin_per_section*round + 1, bin_per_section*(round+1.)); // name, first bin, last bin

  // fit the function and save
  proj_histo[round]->Add(proj_Mtf);
  proj_histo[round]->Fit("gaus", "QEM"); // fit the function to a gaussian

  // Write each fit    
  proj_histo[round]->Write();

  // Get values for Mass Graph
  double real_Mtf_bin_mean = proj_histo[round]->GetMean(); // find mean given in bin graph
  double real_Mtf_bin_std = proj_histo[round]->GetStdDev();
  int num_points = proj_histo[round]->GetEntries(); // gives number of entries as double -> turn to int
  x[round] = (round+0.5)*section_width;  
  y[round] = real_Mtf_bin_mean;
  ex[round] = 0.;
  ey[round] = real_Mtf_bin_std/sqrt(num_points - 1);
  his_P1_Mtf_av->Fill(x[round], y[round]);

  
  std::cout << "P1 Range: [" << section_width*round << ", " << section_width*(round+1.) << "]          "
	    << "no fit std:  " <<  real_Mtf_bin_std << "         no fit mean:  " <<  real_Mtf_bin_mean
	    << "         number of points in sample:  " << num_points << "\n";
  proj_Mtf->Reset("ICESM");
 }
 
// Make mean graph pretty
  his_P1_Mtf_av->SetMarkerStyle(8);  
  his_P1_Mtf_av->SetMarkerSize(1);   
  his_P1_Mtf_av->SetMarkerColor(4);  

 
// Error Bar Graph of Mtf Mean
   TCanvas *c2 = new TCanvas("c2","c2",900,900);// Mass_x_min, Mass_x_max, Mass_y_min, Mass_y_max);
   TGraphAsymmErrors *Mtf_error = new TGraphAsymmErrors(total_sections, x, y, ex, ex, ey, ey);
// Make Pretty and Draw   
   Mtf_error->SetTitle("Mean Field Contribution");            
   Mtf_error->SetLineColor(4);
   Mtf_error->SetMarkerStyle(8);
   Mtf_error->SetMarkerSize(1);                               
   Mtf_error->SetMarkerColor(4);
   Mtf_error->GetYaxis()->SetTitle("Missing Mass [MeV]");
   Mtf_error->GetYaxis()->SetLabelSize(0.03);     
   Mtf_error->GetYaxis()->SetTitleOffset(1.3);                
   Mtf_error->GetYaxis()->CenterTitle(true);                  
   Mtf_error->GetXaxis()->SetTitle("Missing Momentum [MeV]"); 
   Mtf_error->GetXaxis()->SetLabelSize(0.03);     
   Mtf_error->GetXaxis()->CenterTitle(true);                  
   Mtf_error->SetMinimum(860);                                
   Mtf_error->SetMaximum(1040);                               
   Mtf_error->SetFillColor(6);
   Mtf_error->SetFillStyle(3005);
   Mtf_error->Draw("AP");                                     
   c2->Update();
   

// Produce graphs and close
    input_file_mean->Close();
    input_file_SRC->Close();
    output_file->cd();
    his_P1_Mtf->Write();
    his_P1_Mtf_av->Write();
    his_theta_P1prime_q->Write();
    his_Q2_weight->Write();
    Mtf_error->Write();
    
    his_Xb->Write();
    his_Xb_SRC->Write();
    his_Xb_MF ->Write();
    // P1
    his_P1_before->Write();
    his_P1_after->Write();
    his_P1_before_SRC->Write();
    his_P1_after_SRC->Write();
    his_P1_before_MF->Write();
    his_P1_after_MF->Write();
  
    c2->Update();
    output_file->Close();
    return 0;
}


// Add Gaussian Fit:
  /** get gaussian fit
  #include "TF1.h"
  double fit_Mtf_bin_mean = 0;
  double fit_Mtf_bin_std = 0.;
  TF1 *fit = (TF1 *)proj_histo[round]->GetFunction("gaus");
  if(fit != NULL){
    fit_Mtf_bin_mean = fit->GetParameter(1);
    fit_Mtf_bin_std = fit->GetParameter(2);}
  else{
    std:: cout << "No Data to Fit. Set: Mean = STD = 0 \n";} **/
