#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <TStyle.h>

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

//All at the bottom after main. Go there to see what they do.
int which_sector(double phi);
int which_phi(double phi);
int which_theta(double theta);
int which_mom(double mom, double mom_base);
double fitguess(double*x, double*p);
double fiducial_cut_fit(double*x, double*p);
double gaus2(double *x, double *p);

//Everything with _bins can be modified if you want. It changes the range that each histogram covers, so more bins means finer splitting and smaller ranges.
const int mom_bins = 50;
const int sectors = 6;
const int theta_bins = 80;
const int phi_bins = 30;
const double phi_max = 120;
const double theta_max = 160.;

double sector_middle[sectors] = {.001,60.,120.,180.,240.,300.};
double mom_middle[mom_bins];
double phi_middle[sectors][phi_bins];
double theta_middle[theta_bins];

//This is meant to be a program that performs fiducial cuts. It consists of two parts. One that takes in a data file and creates every histogram that we need, and then a separate part (same program) that fits these histograms and outputs the final parameters that represent the cuts. PART 2 TAKES THE OUTPUT OF PART 1. DON'T FORGET THIS.
int main(int argc, char** argv)
{
  gSystem->Load("libTree");
  gStyle->SetOptStat(0);
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Try instead:\n"
           << "\tfiducial run_mode /path/to/output/file /path/to/input/files\n"
           << "for run_mode. 1 = Generate all histograms necessary for fiducial analysis\n"
           << "2 = Perform fiducial analysis on the OUTPUT of run_mode 1.";
      return -1;
    }
  int run_mode = atoi(argv[1]);

  //These arrays store points in each variable. What this does is allow me to loop through events and then create histograms that cover very specific ranges in phi, theta, and momentum.

  double mom_base = .05;
  for (int mom = 0; mom<mom_bins; mom++)
    {
      mom_middle[mom] = mom_base+.1*mom;
    }
  for (int theta = 0; theta<theta_bins; theta++)
    {
      theta_middle[theta] = 1+2*theta;
    }
  for (int sector = 0; sector<sectors; sector++)
    {
      for (int phi =-phi_bins/2; phi<-phi_bins/2+phi_bins; phi++)
        {
          phi_middle[sector][phi+phi_bins/2] = sector_middle[sector]+2*phi+1;
        }
    }

  //This value represents the minimum number of events a histogram must have in order to be analyzed and fit. May get removed, but low statistic histograms are very finicky.
  int threshold = 200;

  //Run mode 1 just makes histograms. It can take a while, but this is why it was separated out. Now I don't need to recreate EVERY histogram in order to fix some bad fits!
  if (run_mode == 1)
    {
      //In case you need to input more than one file, this handles it.
      int numfiles = argc - 3;
      TFile * infile[numfiles];
      cout << "Starting program..." << endl;
      TChain intree("T");
      char * file_name;
      for (int file = 0; file<numfiles; file++)
        {
          file_name = argv[file+3];
          cout << file_name << endl;
          intree.Add(file_name); //Electron and proton Generated
        }

      //This output file will be the INPUT to run mode 2
      TFile * outputFile = new TFile(argv[2], "RECREATE");

      const unsigned int nevents = intree.GetEntries();

      int maxPart = 50;
      int nRun;
      int nParticles;
      /*int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
        double Nu, Q2, Xb, Nu_unc, Q2_unc, Xb_unc, t0;
        double vtx_z_unc [maxPart], vtx_z_cor[maxPart], Mass[maxPart];
        double e_deltat  [maxPart];
        int    stat_sc   [maxPart], stat_ec  [maxPart], stat_dc[maxPart];
        double sc_time   [maxPart], sc_path  [maxPart];
        double ec_time   [maxPart], ec_path  [maxPart];
        double ec_in     [maxPart], ec_out   [maxPart], ec_tot [maxPart];
        double ec_x      [maxPart], ec_y     [maxPart], ec_z   [maxPart];
        double ec_u      [maxPart], ec_v     [maxPart], ec_w   [maxPart];
        double charge    [maxPart], beta     [maxPart];*/
      double mom_x     [maxPart], mom_y    [maxPart], mom_z  [maxPart];
      int Part_type    [maxPart];

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

      intree.SetBranchAddress("nRun"      , &nRun      );
      /*
        intree.SetBranchAddress("nProtons"  , &nProtons  );
        intree.SetBranchAddress("nNeutrons" , &nNeutrons );
        intree.SetBranchAddress("nPiplus"   , &nPiplus   );
        intree.SetBranchAddress("nPiminus"  , &nPiminus  );
        intree.SetBranchAddress("t0"        , &t0        );
        intree.SetBranchAddress("Nu"        , &Nu        );
        intree.SetBranchAddress("Q2"        , &Q2        );
        intree.SetBranchAddress("Xb"        , &Xb        );
        intree.SetBranchAddress("charge"    ,  charge    );
        intree.SetBranchAddress("beta"      ,  beta      );
        intree.SetBranchAddress("vtx_z_unc" ,  vtx_z_unc );
        intree.SetBranchAddress("e_deltat"  ,  e_deltat  );
        intree.SetBranchAddress("stat_sc"   ,  stat_sc   );
        intree.SetBranchAddress("stat_dc"   ,  stat_dc   );
        intree.SetBranchAddress("stat_ec"   ,  stat_ec   );
        intree.SetBranchAddress("sc_time"   ,  sc_time   );
        intree.SetBranchAddress("sc_path"   ,  sc_path   );
        intree.SetBranchAddress("ec_time"   ,  ec_time   );
        intree.SetBranchAddress("ec_path"   ,  ec_path   );
        intree.SetBranchAddress("ec_in"     ,  ec_in     );
        intree.SetBranchAddress("ec_out"    ,  ec_out    );
        intree.SetBranchAddress("ec_tot"    ,  ec_tot    );
        intree.SetBranchAddress("ec_x"      ,  ec_x      );
        intree.SetBranchAddress("ec_y"      ,  ec_y      );
        intree.SetBranchAddress("ec_z"      ,  ec_z      );
        intree.SetBranchAddress("ec_u"      ,  ec_u      );
        intree.SetBranchAddress("ec_v"      ,  ec_v      );
        intree.SetBranchAddress("vtx_z_cor" ,  vtx_z_cor );
        intree.SetBranchAddress("ec_w"      ,  ec_w      );
        intree.SetBranchAddress("Mass"      ,  Mass      );*/
      intree.SetBranchAddress("nParticles", &nParticles);
      intree.SetBranchAddress("mom_x"     ,  mom_x     );
      intree.SetBranchAddress("mom_y"     ,  mom_y     );
      intree.SetBranchAddress("mom_z"     ,  mom_z     );
      intree.SetBranchAddress("Part_type" ,  Part_type );

      // --------------------------------------------------------------------------------------------------
      //Create histograms for each subsection of theta, phi, and momentum based on the central angles defined at the top of the program.
      TH2D * sector_full[sectors];
      TH2D * sector_mom[sectors][mom_bins];
      TH1D * sector_mom_theta[sectors][mom_bins][theta_bins];
      TH1D * sector_mom_phi[sectors][mom_bins][phi_bins];
      char name[100];
      char title[100];
      TH1D * vertex[numfiles];
      cout << "Naming histograms..." << endl;

      //These loops focus on creating and naming the histograms based on the ranges that they cover.
      for(int sector = 0; sector<sectors; sector++)
        {
          double lower = sector_middle[sector]-30.;
          double upper = sector_middle[sector]+30.;
          sprintf(title,"Sector_%d",sector+1);
          sector_full[sector] = new TH2D(title,"",theta_max,0,theta_max,phi_max,lower,upper);

          for(int mom=0; mom<mom_bins; mom++)
            {
              sprintf(title,"Sector_%d_mom_%f_to_%f",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base);
              sector_mom[sector][mom] = new TH2D(title, "",theta_max,0,theta_max,phi_max,lower,upper);

              for(int theta = 0; theta<theta_bins; theta++)
                {
                  int low = theta_middle[theta]-1;
                  int hi = theta_middle[theta]+1;
                  sprintf(title,"Sector_%d_mom_%f_to_%f_theta_%d_to_%d",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base,low,hi);
                  sector_mom_theta[sector][mom][theta] = new TH1D(title, "",phi_max,lower,upper);
                }

              for (int phi = 0; phi<phi_bins; phi++)
                {
                  int low = phi_middle[sector][phi]-1;
                  int hi = phi_middle[sector][phi]+1;
                  sprintf(title,"Sector_%d_mom_%f_to_%f_phi_%d_to_%d",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base,low,hi);
                  sector_mom_phi[sector][mom][phi] = new TH1D(title, "",theta_max,0,theta_max);
                }
            }
        }
      cout << "Finished naming histograms...\n"
           << "Starting loop through events..." << endl;

      double mom_mag = 0;
      double mom_theta = 0;
      double mom_phi = 0;
      TVector3 mom;
      //Loop through events. Meant to now fill the histograms.
      for (int event = 0; event < nevents; event++)
        {
          if (event%100000==0)
            cout << "Working on event " << event << endl;

          intree.GetEvent(event);

          for(int part = 0; part<nParticles; part++)
            {
              //This is where we choose which particle is going to be looked at.
              if ((Part_type[part] != 2212))
                continue;

              TVector3 mom(mom_x[part], mom_y[part], mom_z[part]);
              mom_mag = mom.Mag();
              mom_theta = mom.Theta()*(180/M_PI);
              mom_phi = mom.Phi()*(180/M_PI);
              //Must sanitize phi to stay consistent with the rest of our analyses and those definitions of phi.
              if(mom_phi <-30)
                mom_phi += 360;

              //Check after main for these functions and their explanations. The resultant value will be negative if the function didn't find the right answer. This isn't necessarily bad, it just means that it falls out of our bounds (too high momentum, phi, theta). These regions have very low statistics so they aren't even consider by our fiducial cuts.
              int part_sector = which_sector(mom_phi);
              if (part_sector < 0)
                {
                  cout << "Woah a negative value in which_sector. Check this! It's at event " << event << endl;
                  continue;
                }
              sector_full[part_sector]->Fill(mom_theta,mom_phi);

              int part_mom = which_mom(mom_mag,mom_base);
              if (part_mom<0)
                {
                  cout << "Woah a negative value in which_mom. Check this! It's at event " << event << " mom = " << mom_mag << endl;
                  continue;
                }
              sector_mom[part_sector][part_mom]->Fill(mom_theta,mom_phi);
              int part_theta = which_theta(mom_theta);
              if (part_theta<0)
                {
                  cout << "Woah a negative value in which_theta. Check this! It's at event " << event << " theta = " << mom_theta << endl;
                  continue;
                }
              sector_mom_theta[part_sector][part_mom][part_theta]->Fill(mom_phi);

              int part_phi = which_phi(mom_phi);
              if (part_phi<0)
                {
                  cout << "Woah a negative value in which_phi. Check this! It's at event " << event << endl;
                  continue;
                }
              sector_mom_phi[part_sector][part_mom][part_phi]->Fill(mom_theta);
            }
        }

      cout << "Exited event loop" << endl;
      cout << "Entering write loop" << endl;
      //Just writing out the histograms. The commented out portion is there in case I choose to incorporate phi slices into my analysis. Ignore for now.
      for (int sector = 0; sector<sectors; sector++)
        {
          cout << "In sector " << sector+1 << endl;
          sector_full[sector]->Write();
          delete sector_full[sector];
          //TCanvas * c1[mom_bins];

          for (int mom = 0; mom<mom_bins;mom++)
            {
              cout << "In momentum " << mom_middle[mom] << endl;
              sector_mom[sector][mom]->Write();
              /*char temp[100];
                sprintf(temp, "c%f", mom+1);
                c1[mom] = new TCanvas(temp);
              se[mom]->Print(temp);*/
              delete sector_mom[sector][mom];

              for (int phi = 0; phi<phi_bins;phi++)
                {
                  sector_mom_phi[sector][mom][phi]->Write();
                  delete sector_mom_phi[sector][mom][phi];
                }

              for (int theta = 0; theta<theta_bins;theta++)
                {
                  sector_mom_theta[sector][mom][theta]->Write();
                  delete sector_mom_theta[sector][mom][theta];
                }
            }
        }
      return 0;
    }

  //Run mode 2 goes through the histograms created above, checks to see if they are above threshold and then fits them with our custom fits.
  else if (run_mode == 2)
    {
      //REMEMBER THE INPUT IS THE OUTPUT OF PART 1
      TFile * infile = new TFile(argv[3]);
      cout << "Starting program..." << endl;
      TFile * outfile = new TFile(argv[2], "RECREATE");

      //Get the histograms for each subsection of theta and phi based on the central angles defined at the top of the program as created by run mode 1.
      TH2D * sector_full[sectors];
      TH2D * sector_mom[sectors][mom_bins];
      TH1D * sector_mom_theta[sectors][mom_bins][theta_bins];
      TH1D * sector_mom_phi[sectors][mom_bins][phi_bins];
      cout << "Naming histograms..." << endl;
      char title[256];

      for(int sector = 0;sector<sectors;sector++)
        {
          cout << "In sector " << sector+1 << endl;
          sprintf(title,"Sector_%d",sector+1);
          sector_full[sector] = (TH2D*)infile->Get(title);
          //cout << sector_full[sector] << endl;
          for(int mom=0;mom<mom_bins;mom++)
            {
              sprintf(title,"Sector_%d_mom_%f_to_%f",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base);
              sector_mom[sector][mom] = (TH2D*)infile->Get(title);
              //cout << sector_mom[sector][mom] << endl;
              for(int theta = 0;theta<theta_bins;theta++)
                {
                  int low = theta_middle[theta]-1;
                  int hi = theta_middle[theta]+1;
                  sprintf(title,"Sector_%d_mom_%f_to_%f_theta_%d_to_%d",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base,low,hi);
                  sector_mom_theta[sector][mom][theta] = (TH1D*)infile->Get(title);
                  sector_mom_theta[sector][mom][theta]->GetXaxis()->SetTitle("Phi [degrees]");
                  sector_mom_theta[sector][mom][theta]->GetYaxis()->SetTitle("Counts");

                  sector_mom_theta[sector][mom][theta]->GetXaxis()->CenterTitle();
                  sector_mom_theta[sector][mom][theta]->GetYaxis()->CenterTitle();

                  //cout << sector_mom_theta[sector][mom][theta] << endl;
                }

              for (int phi = 0;phi<phi_bins;phi++)
                {
                  int low = phi_middle[sector][phi]-1;
                  int hi = phi_middle[sector][phi]+1;
                  sprintf(title,"Sector_%d_mom_%f_to_%f_phi_%d_to_%d",sector+1,mom_middle[mom]-mom_base,mom_middle[mom]+mom_base,low,hi);
                  sector_mom_phi[sector][mom][phi] = (TH1D*)infile->Get(title);
                }
            }
        }
      cout << "Finished naming histograms...\n"
           << "Starting fitting process..." << endl;
      double temp[6];

      TF1 * part_boxfit;
      //These arrays store the necessary parameters of our fits for each theta or phi slice.
      double theta_fiducial_upper_bounds[sectors][mom_bins][theta_bins];
      double theta_fiducial_upper_bounds_errors[sectors][mom_bins][theta_bins] = {0};
      double theta_fiducial_lower_bounds[sectors][mom_bins][theta_bins];
      double theta_fiducial_lower_bounds_errors[sectors][mom_bins][theta_bins] = {0};
      //double phi_fiducial_bounds[sectors][mom_bins][phi_bins] {0};
      //double phi_fiducial_bounds_errors[sectors][mom_bins][phi_bins] {0};

      //The arrays above are later fit. These arrays store those secondary parameters. The first array spot of length 2 corresponds to either the upper or lower parts. 0 is upper, 1 is lower.
      double params[2][5][sectors][mom_bins] = {0};
      double params_errors[2][5][sectors][mom_bins] = {0};
      TCanvas *c1[sectors][mom_bins];
      //We choose a sector. We choose a momentum. We go from high to low theta and then look at each slice and fit it with a trapezoidal function. The details are described below.
      for (int sector = 0;sector<sectors;sector++)
        {

          //Everything commented out from here on out is an artifact of potentially looking at phi slices as well, ignore for now.
          /*int theta_bin[mom_bins] = {0};
          int theta_oi[mom_bins] = {0};
          double phi_lower[mom_bins] = {0};
          double phi_upper[mom_bins] = {0};*/

          for (int mom = 0; mom<mom_bins; mom++)
            {
              //theta_bin[mom] = sector_mom[sector][mom]->GetXaxis()->FindBin(30-mom*.6);
              //theta_oi[mom] = theta_bin[mom]/2;
              //cout << theta_bin[mom] << " " << theta_oi[mom] << endl;
              for (int theta = theta_bins-1;theta>=0;theta--)
                {

                  theta_fiducial_upper_bounds[sector][mom][theta] = sector_middle[sector];
                  theta_fiducial_lower_bounds[sector][mom][theta] = sector_middle[sector];

                  if (sector_mom_theta[sector][mom][theta]->Integral(0,theta_max)<threshold)
                    continue;

                  //My attempts at good initial parameters.
                  double max = sector_mom_theta[sector][mom][theta]->GetMaximum();
                  double left_bin = sector_mom_theta[sector][mom][theta]->FindFirstBinAbove(max/5);
                  double right_bin = sector_mom_theta[sector][mom][theta]->FindLastBinAbove(max/5);
                  double left = sector_mom_theta[sector][mom][theta]->GetBinCenter(left_bin);
                  double right = sector_mom_theta[sector][mom][theta]->GetBinCenter(right_bin);

                  //Check below for this function.
                  part_boxfit = new TF1("part_box", fitguess, left-2,right+2,4);

                  part_boxfit->SetParameters(left,1,right,max);
                  //In order to find the proper fit, we limit the parameters to specific ranges. This means trusting our initial parameters.
                  part_boxfit->SetParLimits(0,left-2,left+2);
                  part_boxfit->SetParLimits(2,right-2,right+2);
                  part_boxfit->SetParLimits(1,0,5);
                  sector_mom_theta[sector][mom][theta]->Fit("part_box","q","",left-2,right+2);

                  theta_fiducial_upper_bounds[sector][mom][theta] = part_boxfit->GetParameter(2);
                  //theta_fiducial_upper_bounds_errors[sector][mom][theta] = part_boxfit->GetParameter(1)/2;
                  theta_fiducial_upper_bounds_errors[sector][mom][theta] = 2;
                  theta_fiducial_lower_bounds[sector][mom][theta] = part_boxfit->GetParameter(0);
                  //theta_fiducial_lower_bounds_errors[sector][mom][theta] = part_boxfit->GetParameter(1)/2;
                  theta_fiducial_lower_bounds_errors[sector][mom][theta] = 2;


                  /*if (sector == 5 && mom == 10 && (theta == 32))
                    {
                      theta_fiducial_upper_bounds[sector][mom][theta] = 323.5;
                    }
                  else if (sector == 4 && mom == 7 && (theta == 50))
                    {
                      theta_fiducial_lower_bounds[sector][mom][theta] = 215.5;
                      }*/
                  /*if (sector == 2 && mom == 4 && (theta == 43 || theta == 42))
                    {
                      theta_fiducial_lower_bounds[sector][mom][theta] = 92;
                      theta_fiducial_upper_bounds[sector][mom][theta] = 147;
                    }
                  else if (sector == 3 && mom ==2 && (theta == 53 || theta ==54))
                    theta_fiducial_upper_bounds[sector][mom][theta] = 205;
                  else if (sector == 3 && mom ==2 && (theta == 68 || theta ==69))
                    {
                      theta_fiducial_upper_bounds[sector][mom][theta] = 204.5;
                      theta_fiducial_lower_bounds[sector][mom][theta] = 156.5;
                    }
                  else if (sector == 3 && mom ==2 && (theta == 66 || theta ==69))
                    theta_fiducial_upper_bounds[sector][mom][theta] = 204.5;
                  else if (sector == 3 && mom== 6 && theta == 21)
                  theta_fiducial_lower_bounds[sector][mom][theta] = 155.5;*/
                  /*if (theta==theta_oi[mom])
                    {
                      phi_lower[mom] = sector_mom_theta[sector][mom][theta_oi[mom]]->GetFunction("part_box")->GetParameter(0);
                      phi_upper[mom] = sector_mom_theta[sector][mom][theta_oi[mom]]->GetFunction("part_box")->GetParameter(2);
                    }

                  int phi_bin_lower = sector_mom[sector][mom]->GetYaxis()->FindBin(phi_lower[mom]);
                  int phi_bin_upper = sector_mom[sector][mom]->GetYaxis()->FindBin(phi_upper[mom]);

                  if (phi_bin_lower > 20)
                    {
                      phi_bin_lower = 0;
                      phi_bin_upper = 0;
                    }

                  int phi_oi_lower = phi_bin_lower/2;
                  int phi_oi_upper = phi_bin_upper/2;
                  */

                  sector_mom_theta[sector][mom][theta]->Write();
                  delete sector_mom_theta[sector][mom][theta];
                  /*for (int phi = phi_oi_lower;phi<phi_oi_upper;phi++)
                    {

                      if (sector_mom_phi[sector][mom][phi]->Integral(0,phi_bins*2)<threshold)
                        continue;

                      double sig_left_bin = sector_mom_phi[sector][mom][phi]->FindFirstBinAbove(max/10);
                      double sig_right_bin = sector_mom_phi[sector][mom][phi]->FindLastBinAbove(max/10);
                      double sig_left = sector_mom_phi[sector][mom][phi]->GetBinCenter(sig_left_bin);
                      double sig_right = sector_mom_phi[sector][mom][phi]->GetBinCenter(sig_right_bin);
                      double max_x = sector_mom_phi[sector][mom][phi]->GetBinCenter(sector_mom_phi[sector][mom][phi]->GetMaximumBin());
                      double max = sector_mom_phi[sector][mom][phi]->GetMaximum();

                      TF1 * fit2 = new TF1("fit",gaus2,sig_left-2,sig_right+2,4);
                      fit2->SetParameters(max,max_x,(max_x-sig_left)/2,(sig_right-max_x)/2);
                      //fit2->SetParLimits(2,(max_x-sig_left)/5,(sig_right-max_x)/3);
                      fit2->SetParLimits(0,max/1.2, max*1.2);
                      //fit2->SetParLimits(1,max_x-5,max_x+5);
                      fit2->SetParLimits(2,0,max_x-sig_left);
                      fit2->SetParLimits(3,0,sig_right-max_x);
                      sector_mom_phi[sector][mom][phi]->Fit("fit","q");


                      phi_fiducial_bounds[sector][mom][phi] = fit2->GetParameter(1);
                      phi_fiducial_bounds_errors[sector][mom][phi] = fit2->GetParameter(2);

                      sector_mom_phi[sector][mom][phi]->Write();
                      delete sector_mom_phi[sector][mom][phi];
                      }*/
                }

              int lowtemp = 0;
              int hitemp = 0;
              //lowtemp is the first theta slice (from low to high) that passed the threshold and was fit. Hitemp is the last such slice. These will be part of our initial guesses in the next fit.
              for (int theta = 0; theta<theta_bins;theta++)
                {
                  if ((theta_fiducial_upper_bounds[sector][mom][theta] != sector_middle[sector]) && (lowtemp == 0))
                    {
                      lowtemp = theta-1;
                    }
                }
              for (int theta = theta_bins-1; theta>0;theta--)
                {
                  if ((theta_fiducial_upper_bounds[sector][mom][theta] != sector_middle[sector]) && (hitemp == 0))
                    {
                      hitemp = theta+1;
                      if (hitemp>=55)
                        {
                          hitemp = 55;
                        }
                    }
                }
              double x_error[theta_bins] = {1.};
              char temp[100];
              sprintf(temp, "c%d_%d", sector+1,mom+1);
              c1[sector][mom] = new TCanvas(temp);
              sprintf(temp, "c%d_%d.pdf", sector+1,mom+1);

              //These graphs take all of the boundaries defined in the trapezoids above and puts them together.
              TGraphErrors * theta_lower_gr = new TGraphErrors(theta_bins,theta_middle,theta_fiducial_lower_bounds[sector][mom],x_error,theta_fiducial_lower_bounds_errors[sector][mom]);
              TGraphErrors * theta_upper_gr = new TGraphErrors(theta_bins,theta_middle,theta_fiducial_upper_bounds[sector][mom],x_error,theta_fiducial_upper_bounds_errors[sector][mom]);
              //TGraphErrors * phi_gr = new TGraphErrors(30,phi_fiducial_bounds[sector][mom],(double*)phi_middle[sector],phi_fiducial_bounds_errors[sector][mom],x_error);
              theta_lower_gr->GetXaxis()->SetTitle("Theta [degrees]");
              theta_lower_gr->GetYaxis()->SetTitle("Phi [degrees]");

              theta_lower_gr->GetXaxis()->CenterTitle();
              theta_lower_gr->GetYaxis()->CenterTitle();

              cout << "before hi" << endl;
              //See fit definition after main.
              TF1 * fidfit = new TF1("fidfit",fiducial_cut_fit,0,theta_max, 5);
              TMultiGraph * mg = new TMultiGraph();

              cout << "after hi" << endl;
              mg->Add(theta_upper_gr,"AP");
              mg->Add(theta_lower_gr,"AP");
              //mg->Add(phi_gr,"AP");
              mg->Draw("AP");

              mg->GetXaxis()->SetTitle("Theta [degrees]");
              mg->GetYaxis()->SetTitle("Phi [degrees]");

              mg->GetXaxis()->CenterTitle();
              mg->GetYaxis()->CenterTitle();
              
              //Again, we limit these parameters to help the fit. This however is damaging here because we care about how these parameters evolve over momentum, so we will see...
              /*fidfit->SetParLimits(1,1,10);
              if ((mom == 3) && (sector == 2))
                {
                  fidfit->SetParameters(25,3,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  fidfit->SetParLimits(0,25,25);
                  fidfit->SetParLimits(1,3,3);
                }
              else if((mom==5) && (sector==2))
                {
                  fidfit->SetParameters(25,3,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  fidfit->SetParLimits(1,2.8,3.2);
                }
              else if((mom==3) && (sector==3))
                {
                  fidfit->SetParameters(26,2,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  fidfit->SetParLimits(0,26,26);
                  //fidfit->SetParLimits(1,1.8,2.2);
                }
              else if((mom==4) && (sector==3))
                {
                  fidfit->SetParameters(26,2.5,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  fidfit->SetParLimits(0,26,26);
                }
              else if((mom==2) && (sector==2))
                {
                  fidfit->SetParameters(25,3,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  //fidfit->SetParLimits(1,3,3);
                }
              else if((mom==7) && (sector==1))
                {
                  fidfit->SetParameters(25,4,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  //fidfit->SetParLimits(1,4,4);
                }
              else if((mom==5) && (sector==0))
                {
                  fidfit->SetParameters(25,7,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                  //fidfit->SetParLimits(1,4,4);
                }
              else if((mom==1) && (sector==0))
                {
                  fidfit->SetParameters(25,3,theta_middle[lowtemp]-3,sector_middle[sector],theta_middle[hitemp]);
                  fidfit->SetParLimits(1,2.8,3.2);
                }
              else
                {
                fidfit->SetParameters(25,4,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
                fidfit->SetParLimits(1,1,10);
                }
              */
              fidfit->SetParameters(25,4,theta_middle[lowtemp],sector_middle[sector],theta_middle[hitemp]);
              fidfit->SetParLimits(1,0,10);
                fidfit->SetParLimits(3,sector_middle[sector],sector_middle[sector]);
                //if (sector==0 && mom==1)
                //fidfit->SetParLimits(2,theta_middle[lowtemp]-3.2,theta_middle[lowtemp]-2.8);
                //else
                fidfit->SetParLimits(2,theta_middle[lowtemp]-.2,theta_middle[lowtemp]+.2);

                fidfit->SetParLimits(4,theta_middle[hitemp],theta_middle[hitemp]);
              theta_upper_gr->Fit("fidfit","q","",0,theta_max);

              //As I said params[0] is for the upper bound (higher phi) and [1] is lower (lower phi)
              params[0][0][sector][mom] = fidfit->GetParameter(0);
              params[0][1][sector][mom] = fidfit->GetParameter(1);
              params[0][2][sector][mom] = fidfit->GetParameter(2);
              params[0][3][sector][mom] = fidfit->GetParameter(3);
              params[0][4][sector][mom] = fidfit->GetParameter(4);

              params_errors[0][0][sector][mom] = fidfit->GetParError(0);
              params_errors[0][1][sector][mom] = fidfit->GetParError(1);
              params_errors[0][2][sector][mom] = fidfit->GetParError(2);
              params_errors[0][3][sector][mom] = fidfit->GetParError(3);
              params_errors[0][4][sector][mom] = fidfit->GetParError(4);

              if(params_errors[0][1][sector][mom] > 0)
                {
                  params_errors[0][0][sector][mom] = .2;
                  params_errors[0][1][sector][mom] = .2;
                  params_errors[0][2][sector][mom] = .2;
                  params_errors[0][3][sector][mom] = .2;
                  params_errors[0][4][sector][mom] = .2;
                }
              //cout <<"hi" << sector << " " << mom << endl;
              //The lower curve should mimick the top one fairly closely, so I use the top fit to initialize the bottom.
              /*if (sector ==0 && mom==10)
                {
                  fidfit->SetParameters(-fidfit->GetParameter(0),fidfit->GetParameter(1),fidfit->GetParameter(2),fidfit->GetParameter(3),fidfit->GetParameter(4));
                  //fidfit->SetParLimits(0,-27,-27);
                }
              else if (sector ==2 && mom == 3)
                {
                  fidfit->SetParameters(-26,4,fidfit->GetParameter(2),fidfit->GetParameter(3),fidfit->GetParameter(4));
                  fidfit->SetParLimits(0,-25.8,-26.2);
                  fidfit->SetParLimits(1,3.8,4.2);
                }
              else if (sector ==2 && mom == 5)
                {
                  fidfit->SetParameters(-27,5,fidfit->GetParameter(2),fidfit->GetParameter(3),fidfit->GetParameter(4));
                  fidfit->SetParLimits(0,-27.2,-26.8);
                  fidfit->SetParLimits(1,4.8,5.2);
                }
              else
              {*/
                  fidfit->SetParameters(-fidfit->GetParameter(0),fidfit->GetParameter(1),fidfit->GetParameter(2),fidfit->GetParameter(3),fidfit->GetParameter(4));
                  //}
              theta_lower_gr->Fit("fidfit","q","",0,theta_max);

              //fidfit->SetParameter(1,2);
              //fidfit->SetParLimits(1,2,2);

              params[1][0][sector][mom] = fidfit->GetParameter(0);
              params[1][1][sector][mom] = fidfit->GetParameter(1);
              params[1][2][sector][mom] = fidfit->GetParameter(2);
              params[1][3][sector][mom] = fidfit->GetParameter(3);
              params[1][4][sector][mom] = fidfit->GetParameter(4);

              params_errors[1][0][sector][mom] = fidfit->GetParError(0);
              params_errors[1][1][sector][mom] = fidfit->GetParError(1);
              params_errors[1][2][sector][mom] = fidfit->GetParError(2);
              params_errors[1][3][sector][mom] = fidfit->GetParError(3);
              params_errors[1][4][sector][mom] = fidfit->GetParError(4);

              if(params_errors[1][1][sector][mom] > 0)
                {
                  params_errors[1][0][sector][mom] = .2;
                  params_errors[1][1][sector][mom] = .2;
                  params_errors[1][2][sector][mom] = .2;
                  params_errors[1][3][sector][mom] = .2;
                  params_errors[1][4][sector][mom] = .2;
                }

              theta_upper_gr->SetMarkerStyle(20);
              theta_upper_gr->SetMarkerSize(.85);
              theta_upper_gr->GetFunction("fidfit")->SetLineColor(kBlue);
              theta_lower_gr->SetMarkerStyle(20);
              theta_lower_gr->SetMarkerSize(.85);
              mg->Write();

              //theta_upper_gr->Draw("AP");
              //phi_gr->SetMarkerStyle(20);
              //phi_gr->SetMarkerSize(.85);
              //c1[sector][mom]->Print(temp);
            }

          //These graphs now take the params above and plots and fits them based on momentum dependences.
          double error[mom_bins] = {mom_base};
          double temp_mom[mom_bins];
          for(int i = 0;i<mom_bins;i++)
            {
              temp_mom[i]=mom_middle[i]-mom_base;
            }
          TGraphErrors * a_u_gr = new TGraphErrors(mom_bins,mom_middle,params[0][0][sector],error,params_errors[0][0][sector]);
          TGraphErrors * a_l_gr = new TGraphErrors(mom_bins,mom_middle,params[1][0][sector],error,params_errors[1][0][sector]);
          TGraphErrors * b_u_gr = new TGraphErrors(mom_bins,mom_middle,params[0][1][sector],error,params_errors[0][1][sector]);
          TGraphErrors * b_l_gr = new TGraphErrors(mom_bins,mom_middle,params[1][1][sector],error,params_errors[1][1][sector]);
          TGraphErrors * t_u_gr = new TGraphErrors(mom_bins,mom_middle,params[0][2][sector],error,params_errors[0][2][sector]);
          TGraphErrors * t_l_gr = new TGraphErrors(mom_bins,mom_middle,params[1][2][sector],error,params_errors[1][2][sector]);
          TGraphErrors * s_u_gr = new TGraphErrors(mom_bins,mom_middle,params[0][3][sector],error,params_errors[0][3][sector]);
          TGraphErrors * s_l_gr = new TGraphErrors(mom_bins,mom_middle,params[1][3][sector],error,params_errors[1][3][sector]);
          TGraphErrors * m_u_gr = new TGraphErrors(mom_bins,mom_middle,params[0][4][sector],error,params_errors[0][4][sector]);
          TGraphErrors * m_l_gr = new TGraphErrors(mom_bins,mom_middle,params[1][4][sector],error,params_errors[1][4][sector]);
          a_u_gr->GetXaxis()->SetTitle("Mom");
          a_l_gr->GetXaxis()->SetTitle("Mom");
          b_u_gr->GetXaxis()->SetTitle("Mom");
          b_l_gr->GetXaxis()->SetTitle("Mom");
          t_u_gr->GetXaxis()->SetTitle("Mom");
          t_l_gr->GetXaxis()->SetTitle("Mom");
          s_u_gr->GetXaxis()->SetTitle("Mom");
          s_l_gr->GetXaxis()->SetTitle("Mom");
          m_u_gr->GetXaxis()->SetTitle("Mom");
          m_l_gr->GetXaxis()->SetTitle("Mom");
          /*a_u_gr->SetTitle("Variable a upper");
            a_l_gr->SetTitle("Variable a lower");
            b_u_gr->SetTitle("Variable b upper");
            b_l_gr->SetTitle("Variable b lower");
            t_u_gr->SetTitle("Variable t upper");
            t_l_gr->SetTitle("Variable t lower");
            s_u_gr->SetTitle("Variable s upper");
            s_l_gr->SetTitle("Variable s lower");
            m_u_gr->SetTitle("Variable m upper");
            m_l_gr->SetTitle("Variable m lower");*/
          a_u_gr->SetTitle("");
          a_l_gr->SetTitle("");
          b_u_gr->SetTitle("");
          b_l_gr->SetTitle("");
          t_u_gr->SetTitle("");
          t_l_gr->SetTitle("");
          s_u_gr->SetTitle("");
          s_l_gr->SetTitle("");
          m_u_gr->SetTitle("");
          m_l_gr->SetTitle("");
          TF1 * poly5 = new TF1("poly5","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x",0,1.2);
          a_u_gr->GetXaxis()->SetTitle("Mom [MeV/c] ");
          a_u_gr->GetYaxis()->SetTitle("a");

          a_u_gr->GetXaxis()->CenterTitle();
          a_u_gr->GetYaxis()->CenterTitle();

          a_l_gr->GetXaxis()->SetTitle("Mom [MeV/c]");
          a_l_gr->GetYaxis()->SetTitle("a");

          a_l_gr->GetXaxis()->CenterTitle();
          a_l_gr->GetYaxis()->CenterTitle();

          b_u_gr->GetXaxis()->SetTitle("Mom [MeV/c]");
          b_u_gr->GetYaxis()->SetTitle("b");

          b_u_gr->GetXaxis()->CenterTitle();
          b_u_gr->GetYaxis()->CenterTitle();

          b_l_gr->GetXaxis()->SetTitle("Mom [MeV/c]");
          b_l_gr->GetYaxis()->SetTitle("b");

          b_l_gr->GetXaxis()->CenterTitle();
          b_l_gr->GetYaxis()->CenterTitle();

          t_u_gr->GetXaxis()->SetTitle("Mom [MeV/c]");
          t_u_gr->GetYaxis()->SetTitle("Theta_0");

          t_u_gr->GetXaxis()->CenterTitle();
          t_u_gr->GetYaxis()->CenterTitle();


          a_u_gr->Fit("poly5","q","",.15,.65);
          if (sector == 0)
            {
              temp[0] = a_u_gr->GetFunction("poly5")->GetParameter(0);
              temp[1] = a_u_gr->GetFunction("poly5")->GetParameter(1);
              temp[2] = a_u_gr->GetFunction("poly5")->GetParameter(2);
              temp[3] = a_u_gr->GetFunction("poly5")->GetParameter(3);
              temp[4] = a_u_gr->GetFunction("poly5")->GetParameter(4);
              temp[5] = a_u_gr->GetFunction("poly5")->GetParameter(5);
            }
          if (sector==1 || sector ==2)
            {
              //poly5->SetParameters(temp[0],temp[1],temp[2],temp[3],temp[4],temp[5]);
              //a_u_gr->Fit("poly5","q","",.15,.65);
            }
          cout << a_u_gr->GetFunction("poly5")->GetParameter(0) << "  "
               << a_u_gr->GetFunction("poly5")->GetParameter(1) << "  "
               << a_u_gr->GetFunction("poly5")->GetParameter(2) << "  "
               << a_u_gr->GetFunction("poly5")->GetParameter(3) << "  "
               << a_u_gr->GetFunction("poly5")->GetParameter(4) << "  "
               << a_u_gr->GetFunction("poly5")->GetParameter(5)
               << endl;

          a_l_gr->Fit("poly5","q","",.15,.65);
          cout << a_l_gr->GetFunction("poly5")->GetParameter(0) << "  "
               << a_l_gr->GetFunction("poly5")->GetParameter(1) << "  "
               << a_l_gr->GetFunction("poly5")->GetParameter(2) << "  "
               << a_l_gr->GetFunction("poly5")->GetParameter(3) << "  "
               << a_l_gr->GetFunction("poly5")->GetParameter(4) << "  "
               << a_l_gr->GetFunction("poly5")->GetParameter(5)
               << endl;

          b_u_gr->Fit("poly5","q","",.15,.65);
         
          if (sector ==2)
            {
              //poly5->SetParameters(temp[0],temp[1],temp[2],temp[3],temp[4],temp[5]);
              //poly5->SetParameters(1.5,0,0,0,0,0);
              //b_u_gr->Fit("poly5","q","",.15,.65);
            }
          cout << b_u_gr->GetFunction("poly5")->GetParameter(0) << "  "
               << b_u_gr->GetFunction("poly5")->GetParameter(1) << "  "
               << b_u_gr->GetFunction("poly5")->GetParameter(2) << "  "
               << b_u_gr->GetFunction("poly5")->GetParameter(3) << "  "
               << b_u_gr->GetFunction("poly5")->GetParameter(4) << "  "
               << b_u_gr->GetFunction("poly5")->GetParameter(5)
               << endl;

          b_l_gr->Fit("poly5","q","",.15,.65);
          
          if (sector == 2)
            {
              //              poly5->SetParameters(temp[0],temp[1],temp[2],temp[3],temp[4],temp[5]);
              //poly5->SetParameters(1.5,0,0,0,0,0);
              //b_l_gr->Fit("poly5","q","",.15,.65);
            }
              cout << b_l_gr->GetFunction("poly5")->GetParameter(0) << "  "
               << b_l_gr->GetFunction("poly5")->GetParameter(1) << "  "
               << b_l_gr->GetFunction("poly5")->GetParameter(2) << "  "
               << b_l_gr->GetFunction("poly5")->GetParameter(3) << "  "
               << b_l_gr->GetFunction("poly5")->GetParameter(4) << "  "
               << b_l_gr->GetFunction("poly5")->GetParameter(5)
               << endl;

          t_u_gr->Fit("poly5","q","",.15,.65);
          cout << t_u_gr->GetFunction("poly5")->GetParameter(0) << "  "
               << t_u_gr->GetFunction("poly5")->GetParameter(1) << "  "
               << t_u_gr->GetFunction("poly5")->GetParameter(2) << "  "
               << t_u_gr->GetFunction("poly5")->GetParameter(3) << "  "
               << t_u_gr->GetFunction("poly5")->GetParameter(4) << "  "
               << t_u_gr->GetFunction("poly5")->GetParameter(5)
               << endl;

          a_u_gr->Write();
          a_l_gr->Write();
          b_u_gr->Write();
          b_l_gr->Write();
          t_u_gr->Write();
          t_l_gr->Write();
          s_u_gr->Write();
          s_l_gr->Write();
          m_u_gr->Write();
          m_l_gr->Write();
        }
    

      /*for (int sector=0;sector<sectors;sector++)
        {
          for (int mom=0;mom<mom_bins;mom++)
            {
              for (int theta=0;theta<theta_bins;theta++)
                {
                  outtext << sector << " " << mom << " " << theta_middle[theta] << " ";
                  if (theta_fiducial_lower_bounds[sector][mom][theta] == 0)
                    outtext << "? ? ? ? ";
                  else
                  outtext << theta_fiducial_lower_bounds[sector][mom][theta] << " " << theta_fe2a_ciducial_lower_bounds_errors[sector][mom][theta] << " " << theta_fiducial_upper_bounds[sector][mom][theta] << " " << theta_fiducial_upper_bounds_errors[sector][mom][theta] << " ";

                  outtext << phi_middle[sector][theta] << " ";
                  if ((phi_fiducial_bounds[sector][mom][theta] == 0) || (theta>phi_bins-1))
                    outtext << "? ? " << endl;
                  else
                    outtext << phi_fiducial_bounds[sector][mom][theta] << " " << phi_fiducial_bounds_errors[sector][mom][theta] << endl;
                }
            }
            }*/
      return 0;
    }

  else
    {
      cout << "Wrong run_mode. Please instead follow\n"
           << "1 = Generate all histograms necessary for fiducial analysis\n"
           << "2 = Perform fiducial analysis on the OUTPUT of run_mode 1.";
      return -3;
    }
}

//This function takes in the phi of a specific event. It then loops through the array of central phi and figures out which one it belongs to. Remember that we are looking to place every phi within one degrees of the central one. So a central value of 29 would correspond to phi from 28 to 30.
int which_phi(double phi_value)
{
  for (int phi=0 ; phi<phi_bins ; phi++)
    {
      int part_sector = which_sector(phi_value);
      if (fabs(phi_value-phi_middle[part_sector][phi])<=1)
        {
          return phi;
        }
    }
  return -1;
}

//Same as which_phi, but for theta!
int which_theta(double theta_value)
{
  for (int theta=0 ; theta<theta_bins ; theta++)
    {
      if (fabs(theta_value-theta_middle[theta])<=1)
        {
          return theta;
        }
    }
  return -1;
}

//You get the gist by now.
int which_sector(double phi_value)
{
  for (int sector=0 ; sector<sectors ; sector++)
    {
      if (fabs(phi_value-sector_middle[sector])<=30.)
        {
          return sector;
        }
    }
  return -1;
}

//And definitely by now.
int which_mom(double mom_value, double mom_base)
{
  for (int mom=0;mom<mom_bins;mom++)
    {
      if (fabs(mom_value-mom_middle[mom])<=mom_base)
        {
          return mom;
        }
    }
  return -1;
}

//A double gaussian for the phi slices. Based on Lorenzo Zana's thesis, page 80 and onward.
double gaus2(double *x, double *p)
{
  double A = p[0];
  double m = p[1];
  double s1 = p[2];
  double s2 = p[3];
  if (*x <= m)
    return A*exp(-((*x-m)*(*x-m))/(2*s1*s1));

  if (*x > m)
    return A*exp(-((*x-m)*(*x-m))/(2*s2*s2));

  if (A<=0)
    return 7.E25;

  return 7.E25;
}

//A trapezoid fit.
double fitguess(double*x, double*p)
{
  //Lower bound
  double a = p[0];
  double width = p[1];
  double b = a+width;
  //Upper bound
  double d = p[2];
  double c = d - width;
  double height = p[3];

  if (width < 0)
    return 7.E25;

  if (a>d)
    return 7.E25;

  if ((*x < a) || (*x > d))
    return 0;

  if ((*x > a) && (*x < b))
    return height*(*x-a)/(b-a);

  //if ((*x > b) && (*x < d))
  //return height;
  if ((*x > b) && (*x < c))
    return height;

  if ((*x > c) && (*x < d))
    return height*(d-*x)/(d-c);

  return 7.E25;
}

//Based on Lorenzo Zana's thesis, page 81 on the pdf (that you can google)
double fiducial_cut_fit(double*x, double*p)
{
  double a = p[0];
  double b = p[1];
  double t = p[2];
  double s = p[3];
  double m = p[4];

  //cout << "Sector middle is :" << s << "\n";

  if ((*x < t) || (*x > m))
    return s;
  else
    return s+a*(1.-(1./((*x-t)/b+1.)));

  return s;
}
