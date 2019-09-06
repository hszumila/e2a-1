/*
  This program is a small example to show how a user could analyze
  a skimmed data file and select protons passing leading cuts.
  - Axel 03/8/2018
 */

#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead call:\n"
	   << "\t example_analysis /path/to/skim/file /path/to/output/file\n\n";
      return -1;
    }
  
  // Load the input file
  TFile * infile = new TFile(argv[1]);

  // Set up the tree
  TTree * intree = (TTree*)infile->Get("T");
  const int maxParticles=50;
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles];
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);

    // Produce a histogram, to be filled in the event loop
  TFile * outfile = new TFile(argv[2],"RECREATE");

  TH2D * hist_pmiss_xb = new TH2D("pmiss_xb","Events passing leading cuts;pmiss;Xb;Counts",28,0.3,1.0,40,1,3);

  for (int event =0 ; event < intree->GetEntries() ; event++)
    {
      if (event%10000==0)
	cerr << "Working on event " << event <<"\n";

      intree->GetEvent(event);

      // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      TVector3 q = TVector3(0.,0.,4.461) - mom_e;

      // Loop over particles looking for a proton that will pass leading cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{
	  // Only look at protons
	  if (type[part] != 2212)
	    continue;
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Cut on the angle between p and q
	  if (this_mom.Angle(q) > 25.*M_PI/180.)
	    continue;

	  // Cut on the ratio of p/q
	  double p_over_q = this_mom.Mag() / q.Mag();
	  if (p_over_q < 0.62)
	    continue;
	  if (p_over_q > 0.96)
	    continue;

	  // This proton passes the cuts, it is the leading proton for this event
	  leadingID = part;
	}

      // If we have a proton passing the leading cuts, then fill the histogram
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  TVector3 mom_miss = mom_lead - q;

	  hist_pmiss_xb->Fill(mom_miss.Mag(),xB);
	}
    }

  infile->Close();
  outfile->cd();
  hist_pmiss_xb->Write();
  outfile->Close();

  return 0;
}
