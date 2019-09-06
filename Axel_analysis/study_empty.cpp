#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "constants.h"

using namespace std;

void find_peaks(TH1D *h, vector<double> &list);// Brute force search for peaks, returning a vector of x values
void find_and_fit_all(TH1D * h, vector<double> &posList, vector<double> &errList);

int main(int argc, char ** argv)
{
  if (argc !=3 )
    {
      cerr << "Wrong number of arguments. Instead try\n"
	   << "\tstudy_empty /path/to/input/file /path/to/output/file\n\n";
      return -1;
    }

  TFile * inF = new TFile(argv[1]);
  TTree * inT = (TTree*) inF->Get("data");

  // --------------------------------------------------------------------------------------------------
  // Get the important data from the tree
  const int nEvents = inT->GetEntries();
  const int maxPart = 50;
  int gPart, CCPart, DCPart, ECPart, SCPart, NRun;
  int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
  float STT, W, Yb;
  float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
    SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
    EC_Time[maxPart], EC_Path[maxPart],
    Charge[maxPart], Beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
    pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
    EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], EC_U[maxPart],EC_V[maxPart],EC_W[maxPart], 
    CC_Chi2[maxPart];
  inT->SetBranchAddress("NRun"     ,&NRun   ); // Run number
  inT->SetBranchAddress("gPart"    ,&gPart  ); // Number of particles observed (globally) in the event
  inT->SetBranchAddress("CCPart"   ,&CCPart ); // Number of particles observed in the Cherenkovs
  inT->SetBranchAddress("DCPart"   ,&DCPart ); // Number of particles observed in the Drift Chambers
  inT->SetBranchAddress("ECPart"   ,&ECPart ); // Number of particles observed in the ECal
  inT->SetBranchAddress("SCPart"   ,&SCPart ); // Number of particles observed in the Scintillation Counters (ToFs)
  inT->SetBranchAddress("Stat"     ,Stat    ); // Global status for each particle candidate
  inT->SetBranchAddress("StatDC"   ,StatDC  ); // Drift Chamber status for each particle candidate
  inT->SetBranchAddress("StatCC"   ,StatCC  ); // Cherenkov status for each particle
  inT->SetBranchAddress("StatEC"   ,StatEC  ); // ECal status for each particle
  inT->SetBranchAddress("StatSC"   ,StatSC  ); // Scintillation counter status for each particle  
  inT->SetBranchAddress("particle" ,id_guess); // Guess at the particle ID made by the recon software (maybe not reliable)
  inT->SetBranchAddress("EC_in"    ,EC_in   ); // Inner layer of ECal for each particle
  inT->SetBranchAddress("EC_out"   ,EC_out  ); // Outer layer of ECal for each particle
  inT->SetBranchAddress("EC_tot"   ,EC_tot  ); // Total energy deposit in the ECal for each particle
  inT->SetBranchAddress("Nphe"     ,Nphe    ); // Number of photo-electrons per hit in the Cherenkov detectors
  inT->SetBranchAddress("SC_Time"  ,SC_Time ); // Time in the scintillators per particle
  inT->SetBranchAddress("SC_Path"  ,SC_Path ); // Path Length per particle
  inT->SetBranchAddress("CC_Time"  ,CC_Time ); // Time in the cherenkov per particle
  inT->SetBranchAddress("CC_Path"  ,CC_Path ); // Path Length per particle
  inT->SetBranchAddress("EC_Time"  ,EC_Time ); // Time in the EC per particle
  inT->SetBranchAddress("EC_Path"  ,EC_Path ); // Path Length per particle
  inT->SetBranchAddress("Charge"   ,Charge  ); // Charge per particle
  inT->SetBranchAddress("Beta"     ,Beta    ); // Beta per particle
  inT->SetBranchAddress("Mass"     ,mass    ); // Mass per particle
  inT->SetBranchAddress("Momentum" ,mom     ); // Momentum magnitude per particle
  inT->SetBranchAddress("Momentumx",px      ); // Momentum x component per particle
  inT->SetBranchAddress("Momentumy",py      ); // Momentum y component per particle
  inT->SetBranchAddress("Momentumz",pz      ); // Momentum z component per particle
  inT->SetBranchAddress("Theta"    ,theta   ); // Theta per particle
  inT->SetBranchAddress("Phi"      ,phi     ); // Phi per particle
  inT->SetBranchAddress("TargetZ"  ,targetZ ); // Target Z per particle
  inT->SetBranchAddress("Thetapq"  ,theta_pq); // Angle wrt to q vector per particle
  inT->SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
  inT->SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
  inT->SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
  inT->SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter
  inT->SetBranchAddress("EC_U"     ,EC_U    ); // u positions of hit in the calorimeter
  inT->SetBranchAddress("EC_V"     ,EC_V    ); // v positions of hit in the calorimeter
  inT->SetBranchAddress("EC_W"     ,EC_W    ); // w positions of hit in the calorimeter
  inT->SetBranchAddress("CC_Chi2"  ,CC_Chi2 ); // angle between CC hit and nearest SC hit (in rad)

  // Open the output file and define histograms
  TFile * outF = new TFile(argv[2],"RECREATE");
  TH2D * h2_z_by_sector = new TH2D("z_by_sector","Electrons at 20 deg;z [cm];Sector;Counts",200,-10,10,6,-0.5,5.5);
  TH2D * h2_z_by_phi = new TH2D("z_uncorr","Electrons at 20 deg;z [cm];Phi [deg];Counts",200,-10,10,60,-30.,330.);
  TH2D * h2_z_by_phi_corr = new TH2D("z_corr","Electrons at 20 deg;z [cm];Phi [deg];Counts",200,-10,10,60,-30.,330.);


  for (int event=0 ; event < nEvents ; event++)
    {
      if (event%100000==0)
	cerr << "Working on event " << event << " out of " << nEvents << "\n";

      inT->GetEntry(event);

      if (gPart < 1)
	continue;

      // Require that the electron be between 19 and 21 degrees
      if ((theta[0] < 19.) || (theta[0] > 21.))
	continue;

      int sector = (phi[0]+30.)/60.;
      h2_z_by_sector->Fill(targetZ[0],sector);
      h2_z_by_phi->Fill(targetZ[0], phi[0]);

    }

  // Loop over sector and do fits
  double us_foil_z[6];
  double us_foil_dz[6];
  double phis[6]={0.,M_PI/3.,2.*M_PI/3.,M_PI,4.*M_PI/3.,5.*M_PI/3.};
  double dphis[6]={0.1,0.1,0.1,0.1,0.1,0.1};

  // Fit the upstream foil position (cleanest)
  TF1 us_foil("us_foil","gaus",-10,-5);
  for (int s=0 ; s<6 ; s++)
    {
      char temp[100];
      sprintf(temp,"z_%d",s);
      TH1D * h1_z = h2_z_by_sector->ProjectionX(temp,s+1,s+1);
     
      vector<double> peaks;
      vector<double> errs;
      find_and_fit_all(h1_z,peaks,errs);

      if (peaks.size() > 0)
	{
	  us_foil_z[s] = peaks[0];
	  us_foil_dz[s] = errs[0];
	}

      h1_z->Write();
    }

  // Form a t-graph errors
  TGraphErrors * us_foil_fits = new TGraphErrors(6,phis,us_foil_z,dphis,us_foil_dz);
  TF1 sinfit("sinfit","[0]*sin(x - [1]) + [2]",0.,2.*M_PI);
  us_foil_fits->Fit(&sinfit,"Q");
  us_foil_fits->Write();

  double amp = sinfit.GetParameter(0);
  double phase = sinfit.GetParameter(1);
 
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event%100000==0)
	cerr << "Working on event " << event << " out of " << nEvents << "\n";
      
      inT->GetEntry(event);
      
      if (gPart < 1)
	continue;
      
      // Require that the electron be between 19 and 21 degrees
      if ((theta[0] < 19.) || (theta[0] > 21.))
	continue;
      
      double z_corr = targetZ[0] - amp*sin(phi[0]*M_PI/180. - phase);
      
      h2_z_by_phi_corr->Fill(z_corr, phi[0]);
    }

  
  TH1D * h1_zcorr = h2_z_by_phi_corr->ProjectionX("z_corr_1d",1,60);
  vector<double> peaks;
  vector<double> errs;
  find_and_fit_all(h1_zcorr,peaks,errs);
  
  cout << amp << " " << phase*180./M_PI;
  for (unsigned int i=0 ; i < peaks.size() ; i++)
    {
      cout << " " << peaks[i] << " " << errs[i];
    }
  cout << "\n";

  h1_zcorr->Write();

  inF->Close();
  h2_z_by_sector->Write();
  h2_z_by_phi->Write();
  h2_z_by_phi_corr->Write();
  outF->Close();
  
  return 0;
}

// Up-stream and downstream window fit functions are not currently used.
double us_window(double *x, double *p)
{
  double a=p[0];
  double m=p[1];
  double s=p[2];
  double b=p[3];

  double temp=(*x - m)/s;
  return a*exp(-0.5*temp*temp) + b/(1. + exp(-(*x-m)/0.5));
}

double ds_window(double *x, double *p)
{
  double a=p[0];
  double m=p[1];
  double s=p[2];
  double b=p[3];

  double temp=(*x - m)/s;
  return a*exp(-0.5*temp*temp) + b/(1. + exp(*x/0.5));
}

// Brute force search for peaks, returning a vector of bin values
void find_peaks(TH1D *h, vector<double> &list)
{
  list.clear();
  int nBins = h->GetNbinsX();
  double maxVal = h->GetMaximum();
  for (int i=1 ; i<= nBins ; i++)
    {
      double thisVal = h->GetBinContent(i);
      bool biggerThanLNeighbors=true;
      bool biggerThanRNeighbors=true;
      for (int j=1 ; j<=5 ; j++)
	{
	  if (i-j>0)
	    biggerThanLNeighbors &= (thisVal >= h->GetBinContent(i-j));
	  if (i+j<=nBins)
	    biggerThanRNeighbors &= (thisVal >= h->GetBinContent(i+j));
	}

      if (biggerThanLNeighbors && biggerThanRNeighbors && (thisVal > 0.25*maxVal))
	list.push_back(h->GetBinCenter(i));
    }
}

void find_and_fit_all(TH1D * h, vector<double> &posList, vector<double> &errList)
{
  vector<double> peakList;
  find_peaks(h,peakList);
  
  for (unsigned int i=0 ; i < peakList.size() ; i++)
    {
      TFitResultPtr r = h->Fit("gaus","Q+S","",peakList[i]-0.7,peakList[i]+0.7);
      posList.push_back(r->Parameter(1));
      errList.push_back(r->Error(1));
    }
}
