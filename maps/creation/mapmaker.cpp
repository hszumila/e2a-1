#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>
#include "TSystem.h"

#include "Fiducial.h"
#include "constants.h"
#include "global_variables.h"


using namespace std;

const int pbins = 100;
const int costbins = 200;
const int phibins = 360;

TH3D * generated;
TH3D * accepted;

int main(int argc, char ** argv)
{
  if (argc < 7)
    {
      cerr << "Wrong number of arguments. Instead try\n"
           << "\tmapmaker /path/to/inputs /path/to/output energy torus_current target particle_of_interest\n\n";
      cerr << "For particle of interest, please use either e, p, pip, or pim for electron, proton, pi-plus, or pi-minus respectively" << endl;
      cerr << "For energy use either 1161, 2261, or 4461" << endl;
      cerr << "For target use either 3, 4, or 12 for He3, He4 or, C12" << endl;
      cerr << "For torus_current use either 2250 for energies of 2261 and 4461. Use 1500 or 750 for energy of 1161." << endl;
      return -3;
    }
  int numfiles = argc-6;
  TFile * infile[numfiles];
  for (int i = 0; i < numfiles; i++)
    {
      infile[i] = new TFile(argv[i+1]);
    }
  string tab_targ;
  double tab_E1 = atoi(argv[numfiles+2]);
  double tab_torus = atoi(argv[numfiles+3]);
  double tab_mini = 5996;

  if (atoi(argv[numfiles+4]) == 3)
    tab_targ = "3He";
  else  if (atoi(argv[numfiles+4]) == 12)
    tab_targ = "12C";
  else if (atoi(argv[numfiles+4]) == 4)
    tab_targ = "4He";
  else
    {
      cerr << "Wrong target type, choose either 3, 4, or 12 for 3He, 4He, or 12C" << endl;
      return -4;
    }

  string particle_oi = argv[numfiles+5];

  if (particle_oi != "e" && particle_oi != "p" && particle_oi != "pip" && particle_oi != "pim")
    {
      cout << "Wrong particle of interest, please use either e, p, pip, or pim for electron, proton, pi-plus, or pi-minus respectively" << endl;
      return -5;
    }

  cout << "Before fiducial class" << endl;
  Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, true);  // Create an instance of the Fiducial Class
  const double EC_in_cut = fid_params.EC_in_cut();
  const double el_EC_cut = fid_params.el_EC_cut();

  cout << "After fiducial class" << endl;

  cout << tab_E1 << " " << tab_torus << " " << tab_targ << endl;
  cout << atoi(argv[numfiles+4]) << endl;
  TFile * outfile = new TFile(argv[numfiles+1],"RECREATE");

  generated = new TH3D("Generated Particles","Generated Particles",pbins,0,5,costbins,-1,1,phibins,-30,330);
  accepted = new TH3D("Accepted Particles","Accepted Particles",pbins,0,5,costbins,-1,1,phibins,-30,330);


  const int maxPart = 50;
  int gPart, CCPart, DCPart, ECPart, SCPart, NRun;
  int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
  float STT, W, Yb;
  float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
    SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
    EC_Time[maxPart], EC_Path[maxPart],
    charge[maxPart], beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
    pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
    EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], EC_U[maxPart],EC_V[maxPart],EC_W[maxPart], 
    CC_Chi2[maxPart];
  float targetZ_g[maxPart], theta_g[maxPart], phi_g[maxPart], mom_g[maxPart], px_g[maxPart], py_g[maxPart], pz_g[maxPart];
  int particle_g[maxPart];
  int num_g;
  for (int file = 0; file < numfiles; file++)
    {
      TTree * intree = (TTree*)infile[file]->Get("data");
      if (!intree)
        {
          cerr << "Tree not found!!!\n";
          return -3;
        }
      intree->SetBranchAddress("gPart",&gPart);
      intree->SetBranchAddress("StatCC",StatCC);
      intree->SetBranchAddress("StatDC",StatDC);
      intree->SetBranchAddress("StatEC",StatEC);
      intree->SetBranchAddress("StatSC",StatSC);
      intree->SetBranchAddress("Stat",Stat);
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
      intree->SetBranchAddress("Number_g",&num_g);
      intree->SetBranchAddress("particle_g",particle_g);
      intree->SetBranchAddress("TargetZ_g",targetZ_g);
      intree->SetBranchAddress("Theta_g",theta_g);
      intree->SetBranchAddress("Phi_g",phi_g);
      intree->SetBranchAddress("Momentum_g",mom_g);
      intree->SetBranchAddress("Momentumx_g",px_g);
      intree->SetBranchAddress("Momentumy_g",py_g);
      intree->SetBranchAddress("Momentumz_g",pz_g);

      for(int event = 0; event<intree->GetEntries();event++)
        {
          intree->GetEvent(event);
          if (event%100000==0)
            cout << "File " << file+1 << " and event " << event << " out of " << intree->GetEntries() << endl;
	      
	  if (particle_oi == "e")
	    {
	      if (num_g <= 0)
		continue;

	      double cost_g = TMath::Cos(theta_g[0]*M_PI/180);
	      int sec_g = (phi_g[0]+30.)/60;

	      generated->Fill(mom_g[0],cost_g,phi_g[0]);

	      if (gPart < 1)
		continue;

	      if (!(			(StatEC[0] > 0) && // EC status is good for the electron candidate
					(StatDC[0] > 0) && // DC status is good for the electron candidate
					(StatCC[0] > 0) && // CC status is good for the electron candidate
					(StatSC[0] > 0) && // SC status is good for the electron candidate
					(charge[0] < 0)    // Electron candidate curvature direction is negative
					))
		{continue;}

	      double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC

	      // Electron E/p cut
	      // We no longer due this cut due to the simulation mismodeling the 
	      // calorimeter sampling fraction. 
	      // if (!( (EC_in [0] > EC_in_cut) &&  // Electron candidate has enough energy deposit in inner layer of EC
	      //        (el_cand_EC > el_EC_cut)  &&  // Enough total energy in the EC
	      //        (fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)))) // Electron PID (E/p)
	      // {continue;}

	      // Electron fiducial cuts
	      // We no longer do this cut due to problems with the granularity of
	      // the map binning.
	      TVector3 T3_e_mom(px[0],py[0],pz[0]);
	      TVector3 e_ec_xyz(EC_X[0],EC_Y[0],EC_Z[0]);
	      // Electron Fiducial cuts
	      //if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
	      //if (!fid_params.CutUVW_e(e_ec_xyz)       ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);
	      TVector3 T3_e_mom_g(px_g[0],py_g[0],pz_g[0]);
	      //if (!fid_params.e_inFidRegion(T3_e_mom_g)) continue; // Electron theta-phi cut
	      
	      // Require that the generated and reconstructed sectors match.
	      int sec = (phi[0]+30.)/60;
	      if (sec != sec_g)
		continue;

	      accepted->Fill(mom_g[0],cost_g,phi_g[0]);
	    }
	  else if (particle_oi == "p" || particle_oi == "pip")
	    {
	      if (num_g < 2)
		continue;

	      double cost_g = TMath::Cos(theta_g[1]*M_PI/180);

	      if (!((StatEC[0] > 0) && // EC status is good for the electron candidate
		    (StatDC[0] > 0) && // DC status is good for the electron candidate
		    (StatCC[0] > 0) && // CC status is good for the electron candidate
		    (StatSC[0] > 0) && // SC status is good for the electron candidate
		    (charge[0] < 0)    // Electron candidate curvature direction is negative
		    ))
		{continue;}

	      // We have a good electron candidate. Go ahead and fill generated
	      generated->Fill(mom_g[1],cost_g,phi_g[1]);

	      // Define the sector
	      int sec_g = (phi_g[1]+30.)/60;
	      if (gPart <= 1)
		continue;

	      bool pass = false;

	      for (int part = 1; part < gPart; part++)
		{
		  //Positive particle test
		  if (!(
			(StatSC[part] > 0) && 		// SC status is good for the positive candidate
			(StatDC[part] > 0) &&              // DC status is good for the positive candidate
			(Stat  [part] > 0) &&		// Global status is good for the positive candidate
			(charge[part] > 0) 		// Charge is positive
			))
		    {continue;}

		  TVector3 T3_p_mom(px[part],py[part],pz[part]);
		  double e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;
		  
		  if (particle_oi == "p")
		    {
		      double beta_assuming_proton = mom[part]/sqrt(mom[part]*mom[part] + mP*mP);
		      double p_t0 = SC_Time[part] - SC_Path[part]/(beta_assuming_proton * c_cm_ns);
		      double delta_t = p_t0 - e_t0;
		  
		      //if (!fid_params.in_p_deltaT(delta_t, mom[part], pdeltat_sig_cutrange)) {continue;} // Proton PID (delta T vs p)
		      //if (!fid_params.pFiducialCut(T3_p_mom)) {continue;}

		      TVector3 T3_p_mom_g(px_g[1],py_g[1],pz_g[1]);
		      //if (!fid_params.pFiducialCut(T3_p_mom_g)) {continue;}
		      int sec = (phi[part]+30.)/60;

		      if (sec != sec_g)
			continue;

		      pass = true;
		    }		  
		  else if (particle_oi == "pip")
		    {
		      double beta_assuming_pip = mom[part]/sqrt(mom[part]*mom[part] + mpc*mpc);
		      double pip_t0 = SC_Time[part] - SC_Path[part]/(beta_assuming_pip * c_cm_ns);
		      double pip_delta_t = pip_t0 - e_t0;

		      //if (!fid_params.in_pip_deltaT(pip_delta_t, mom[part], pipdeltat_sig_cutrange)) {continue;} // Proton PID (delta T vs p)
		      //if (!fid_params.pFiducialCut(T3_p_mom)) {continue;}

		      int sec = (phi[part]+30.)/60;

		      if (sec != sec_g)
			continue;

		      pass = true;
		    }
		}
	      if (pass)
		accepted->Fill(mom_g[1],cost_g,phi_g[1]);
	    }	    
	}
    }
  for (int p = 1; p<=pbins;p++)
    {
      for (int phi = 1; phi<=phibins; phi++)
        {
          for (int cost = 1; cost<=costbins; cost++)
            {
              if (generated->GetBinContent(p,cost,phi)==0)
                generated->SetBinContent(p, cost, phi, 10);
            }
        }
    }
  generated->Write();
  accepted->Write();
  outfile->Close();

  return 0;
}
