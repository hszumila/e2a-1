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

using namespace std;

const int pbins = 100;
const int costbins = 200;
const int phibins = 360;

TH3D * generated;
TH3D * accepted;

int main(int argc, char ** argv)
{
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Instead try\n"
           << "\tmapmaker /path/to/inputs /path/to/output particle\n\n";
      cerr << "For particle of interest, please use e, p, pip, or pim, for electron, proton, pi+, or pi- respectively.\n\n";
      return -1;
    }
  int numfiles = argc-3;

  TFile * outfile = new TFile(argv[numfiles+1],"RECREATE");
  string particle_oi = argv[numfiles+2];
  string targ("deut");

  if ((particle_oi != "e") && (particle_oi != "p") && (particle_oi != "pip") && (particle_oi != "pim"))
    {
      cout << "Wrong particle of interest, please use either e, p, pip, or pim for electron, proton, pi-plus, or pi-minus respectively" << endl;
      return -5;
    }
  else
    {
      char temp_name[100];
      char temp_title[100];

      sprintf(temp_name,"%s_%s_gen",targ.c_str(),particle_oi.c_str());
      sprintf(temp_title,"Generated %s;momentum [GeV/c];Cos(theta);Phi [deg]",particle_oi.c_str());
      generated = new TH3D(temp_name,temp_title,pbins,0,5,costbins,-1,1,phibins,-30,330);
      
      sprintf(temp_name,"%s_%s_acc",targ.c_str(),particle_oi.c_str());
      sprintf(temp_title,"Accepted %s;momentum [GeV/c];Cos(theta);Phi [deg]",particle_oi.c_str());
      accepted = new TH3D(temp_name,temp_title,pbins,0,5,costbins,-1,1,phibins,-30,330);
    }

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
  int num_g, particle_g[maxPart], particle[maxPart];

  for (int file = 0; file < numfiles; file++)
    {
      TFile * thisFile = new TFile(argv[file+1]);

      TTree * intree = (TTree*)thisFile->Get("data");
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
      intree->SetBranchAddress("particle",particle);
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
	      if (num_g < 1)
		{
		  cerr << "Error: event has no generated particles!\n";
		  continue;
		}
	      
	      // Fill the generated histogram based on generated values
              double cost_g = TMath::Cos(theta_g[0]*M_PI/180);
              generated->Fill(mom_g[0],cost_g,phi_g[0]);
	      
	      // Calculate the generated electron sector, crucial for later.
              int sec_g = (phi_g[0]+30.)/60;
	      
	      // If there are no reconstructed particles, that's fine, no acceptance
	      if (gPart < 1)
		continue;
	      
	      // Assume that index [0] is the electron candidate
              if (!((StatEC[0] > 0) && // EC status is good for the electron candidate
		    (StatDC[0] > 0) && // DC status is good for the electron candidate
		    (StatCC[0] > 0) && // CC status is good for the electron candidate
		    (StatSC[0] > 0) && // SC status is good for the electron candidate
		    (charge[0] < 0)    // Electron candidate curvature direction is negative
		    ))
                {continue;}

	      // Compare to the generated sector
              int sec = (phi[0]+30.)/60;
	      
              if (sec != sec_g)
                continue;
	      
	      // Particle met acceptance criteria
              accepted->Fill(mom_g[0],cost_g,phi_g[0]);
            }
          else if(particle_oi =="p")
            {
	      // Determine if the electron is good
              if (gPart < 1)
                continue;
	      
              if (!( (StatEC[0] > 0) && // EC status is good for the electron candidate
                     (StatDC[0] > 0) && // DC status is good for the electron candidate
                     (StatCC[0] > 0) && // CC status is good for the electron candidate
                     (StatSC[0] > 0) && // SC status is good for the electron candidate
                     (charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
                continue;
	      
	      // The electron is good!
	      
	      if (num_g < 2)
		{
		  cerr << "Error: event does not have two generated particles!\n";
		  continue;
		}
	      
	      // Define the relevant generated proton quantities, and fill the generated histogram
              double cost_g = TMath::Cos(theta_g[1]*M_PI/180);
              generated->Fill(mom_g[1],cost_g,phi_g[1]);
	      
              // Define a sector variable
              int sec_g = (phi_g[1]+30.)/60;
	      
	      // Do we have a proton passing selection criteria?	      
              bool pass = false;
              for (int part = 1; (part < gPart)&&(!pass); part++)
                {
		  //cerr << "Working on particle " << part << " out of " << gPart << "\n";
		  
                  //Positive particle test
                  if (!(
                        (StatSC[part] > 0) && 		// SC status is good for the positive candidate
                        (StatDC[part] > 0) &&              // DC status is good for the positive candidate
                        (Stat  [part] > 0) &&		// Global status is good for the positive candidate
                        (charge[part] > 0) 		// Charge is positive
                        ))
                    {continue;}
		  
                  int sec = (phi[part]+30.)/60;
		  
                  if (sec != sec_g)
		    {
		      //cerr << "Generated in sector " << sec_g << " but reconstructed in sector " << sec << "\n";
		      continue;
		    }
		  
                  //if (particle[part] != 2212)
		  // {
		  //cerr << "Particle guess is " << particle[part] << "\n";
		  //    continue;
		  //  }
		  
		  // Looks like a suitable proton
                  pass = true;
                }
              if (pass)
                accepted->Fill(mom_g[1],cost_g,phi_g[1]);
            }
	  // ************************************************************************
	  // *********** PI+ ********************************************************
	  // ************************************************************************
	  else if(particle_oi =="pip") 
	    {
	      // Determine if the electron is good
	      if (gPart < 1)
                continue;
	      
	      if (!( (StatEC[0] > 0) && // EC status is good for the electron candidate
		     (StatDC[0] > 0) && // DC status is good for the electron candidate
		     (StatCC[0] > 0) && // CC status is good for the electron candidate
		     (StatSC[0] > 0) && // SC status is good for the electron candidate
		     (charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
                continue;
	      
	      // The electron is good!
	      
	      if (num_g < 2)
		{
		  cerr << "Error: event does not have two generated particles!\n";
		  continue;
		}

	      // Define the relevant generated proton quantities, and fill the generated histogram
              double cost_g = TMath::Cos(theta_g[1]*M_PI/180);
              generated->Fill(mom_g[1],cost_g,phi_g[1]);
	      
              // Define a sector variable
              int sec_g = (phi_g[1]+30.)/60;

	      // Do we have a pi+ passing selection criteria?	      
              bool pass = false;
              for (int part = 1; (part < gPart)&&(!pass); part++)
                {
		  //cerr << "Working on particle " << part << " out of " << gPart << "\n";

                  //Positive particle test
                  if (!(
                        (StatSC[part] > 0) && 		// SC status is good for the positive candidate
                        (StatDC[part] > 0) &&           // DC status is good for the positive candidate
                        (Stat  [part] > 0) &&		// Global status is good for the positive candidate
                        (charge[part] > 0) 		// Charge is positive
                        ))
                    {continue;}

                  int sec = (phi[part]+30.)/60;
		  
                  if (sec != sec_g)
		    {
		      //cerr << "Generated in sector " << sec_g << " but reconstructed in sector " << sec << "\n";
		      continue;
		    }
		  
		  // Looks like a suitable pi+
                  pass = true;
		  break;
                }
	      
              if (pass)
                accepted->Fill(mom_g[1],cost_g,phi_g[1]);
            }
	  // ************************************************************************
	  // *********** PI- ********************************************************
	  // ************************************************************************
	  else if(particle_oi =="pim")
	    {
	      // Determine if the electron is good
	      if (gPart < 1)
                continue;
	      
	      if (!( (StatEC[0] > 0) && // EC status is good for the electron candidate
		     (StatDC[0] > 0) && // DC status is good for the electron candidate
		     (StatCC[0] > 0) && // CC status is good for the electron candidate
		     (StatSC[0] > 0) && // SC status is good for the electron candidate
		     (charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
                continue;
	      
	      // The electron is good!
	      
	      if (num_g < 2)
		{
		  cerr << "Error: event does not have two generated particles!\n";
		  continue;
		}

	      // Define the relevant generated proton quantities, and fill the generated histogram
              double cost_g = TMath::Cos(theta_g[1]*M_PI/180);
              generated->Fill(mom_g[1],cost_g,phi_g[1]);
	      
              // Define a sector variable
              int sec_g = (phi_g[1]+30.)/60;

	      // Do we have a pi- passing selection criteria?	      
              bool pass = false;
              for (int part = 1; (part < gPart)&&(!pass); part++)
                {
		  //cerr << "Working on particle " << part << " out of " << gPart << "\n";

                  //Negative particle test
                  if (!(
                        (StatSC[part] > 0) && 		// SC status is good for the candidate
                        (StatDC[part] > 0) &&           // DC status is good for the candidate
                        (Stat  [part] > 0) &&		// Global status is good for the candidate
                        (charge[part] < 0) 		// Charge is negative
                        ))
                    {continue;}

                  int sec = (phi[part]+30.)/60;
		  
                  if (sec != sec_g)
		    {
		      //cerr << "Generated in sector " << sec_g << " but reconstructed in sector " << sec << "\n";
		      continue;
		    }
		  
		  // Looks like a suitable pi-
                  pass = true;
		  break;
                }
	      
              if (pass)
                accepted->Fill(mom_g[1],cost_g,phi_g[1]);
            }
        }
      thisFile->Close();
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
  outfile->cd();
  generated->Write();
  accepted->Write();

  outfile->Close();

  return 0;
}
