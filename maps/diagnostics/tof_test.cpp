#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>
#include "TSystem.h"

#include "Fiducial.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

int main(int argc, char ** argv)
{
  gSystem->Load("libTree");
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Try instead:\n"
           << "./tof_test /path/to/output /path/to/input\n\n";
      return -1;
    }

  double tab_E1 = 2261;
  double tab_torus = 2250;
  double tab_mini = 5996;
  string tab_targ = "12C";

  Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, true);
  TFile * infile = new TFile(argv[2]);
  TFile * outfile = new TFile(argv[1], "RECREATE");

  TTree * intree = (TTree*)infile->Get("data");
  const int nevents = intree->GetEntries();

  const int maxPart = 50;
	int gPart, CCPart, DCPart, ECPart, SCPart, NRun;
	int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], SC_Pad[maxPart], id_guess[maxPart];
	float STT, W, Yb;
	float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
    SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
    EC_Time[maxPart], EC_Path[maxPart],
    charge[maxPart], beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
    pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
    EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], EC_U[maxPart],EC_V[maxPart],EC_W[maxPart],
    CC_Chi2[maxPart];

      intree->SetBranchAddress("gPart",&gPart);
      intree->SetBranchAddress("StatCC",StatCC);
      intree->SetBranchAddress("StatDC",StatDC);
      intree->SetBranchAddress("StatEC",StatEC);
      intree->SetBranchAddress("StatSC",StatSC);
      intree->SetBranchAddress("Charge",charge);
      intree->SetBranchAddress("EC_in",EC_in);
      intree->SetBranchAddress("EC_out",EC_out);
      intree->SetBranchAddress("EC_tot",EC_tot);
      intree->SetBranchAddress("EC_U",EC_U);
      intree->SetBranchAddress("EC_V",EC_V);
      intree->SetBranchAddress("EC_W",EC_W);
      intree->SetBranchAddress("EC_X",EC_X);
      intree->SetBranchAddress("EC_Y",EC_Y);
      intree->SetBranchAddress("EC_Z",EC_Z);
      intree->SetBranchAddress("Nphe",Nphe);
      intree->SetBranchAddress("Mass",mass);
      intree->SetBranchAddress("Beta",beta);
      intree->SetBranchAddress("TargetZ",targetZ);
      intree->SetBranchAddress("Theta",theta);
      intree->SetBranchAddress("Phi",phi);
      intree->SetBranchAddress("Momentum",mom);
      intree->SetBranchAddress("Momentumx",px);
      intree->SetBranchAddress("Momentumy",py);
      intree->SetBranchAddress("Momentumz",pz);
      intree->SetBranchAddress("SC_Time",SC_Time);
      intree->SetBranchAddress("SC_Path",SC_Path);
      intree->SetBranchAddress("SC_Pad",SC_Pad);

  TH2D * tof_phi = new TH2D("tof_vs_phi", "tof_vs_phi",660,0,66000,360,-30,330);
  TH2D * tof_theta = new TH2D("tof_vs_theta", "tof_vs_theta",660,0,66000,180,0,180);
  TH1D * tof = new TH1D("tof", "tof",660,0,660);
  for (int event = 0;event<nevents;event++)
    {
      if (event%100000==0)
        cout << "Event " << event << " of " << nevents << endl;
      intree->GetEvent(event);
      for (int part = 1; part < gPart; part++)
        {
          //Positive particle test
          if (!(           StatSC[part] > 0) && 		// SC status is good for the positive candidate
              (StatDC[part] > 0) &&              // DC status is good for the positive candidate
              (Stat  [part] > 0) &&		// Global status is good for the positive candidate
              (charge[part] > 0) 		// Charge is positive
              )
            {continue;}

          TVector3 T3_p_mom(px[part],py[part],pz[part]);
          double e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;
          double beta_assuming_proton = mom[part]/sqrt(mom[part]*mom[part] + mP*mP);
          double p_t0 = SC_Time[part] - SC_Path[part]/(beta_assuming_proton * c_cm_ns);
          double delta_t = p_t0 - e_t0;

          if (!fid_params.in_p_deltaT(delta_t, mom[part], pdeltat_sig_cutrange)) {continue;} // Proton PID (delta T vs p)
          if (!fid_params.pFiducialCut(T3_p_mom)) {continue;}
          tof_phi->Fill(SC_Pad[part],phi[part]);
          tof_theta->Fill(SC_Pad[part],theta[part]);
          int tof_num = SC_Pad[part]/100;
          tof->Fill(tof_num);
        }
    }
  tof->Write();
  delete tof;
  tof_phi->Write();
  delete tof_phi;
  tof_theta->Write();
  delete tof_theta;
  outfile->Write();
  delete outfile;
  return 0;
}
