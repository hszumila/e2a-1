// ======================================================
// CODE TO DETERMINE PARAMETERS NEEDED FOR ELECTRON PID
// CUTS IN THE CASE OF SIMULATIONS
// ======================================================
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <TStyle.h>
#include"TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>

#include <TGraphErrors.h>
#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\tsim_pid_electron /path/to/input/file Ebeam_MeV nuclear_target\n\n";
		return -1;	
	}

	// --------------------------------------------------------------------------------------------------
	// Open up the Root file
	TFile * f = new TFile(argv[1]);
	if (f)
		cerr << "Successfully opened file " << argv[1] << "\n";
	else
	{
		cerr << "Could not open file " << argv[1] << "\n\tExiting...\n";
		return -2;
	}

	istringstream iss (argv[2]);
	int in_num_part;
	iss >> in_num_part;

	// --------------------------------------------------------------------------------------------------
	// Open up the tree, and get the important data
	TTree * t = (TTree*)f->Get("data");
	const int nEvents = t->GetEntries();
	const int maxPart = 50;
	int gPart, CCPart, DCPart, ECPart, SCPart;
	int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
	float STT, W, Yb;
	float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
	      SC_Time[maxPart], SC_Path[maxPart],
	      Charge[maxPart], Beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
	      pz[maxPart], theta[maxPart], phi[maxPart], theta_pq[maxPart],
	      EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart];

	t->SetBranchAddress("gPart"    ,&gPart  ); // Number of particles observed (globally) in the event
	t->SetBranchAddress("CCPart"   ,&CCPart ); // Number of particles observed in the Cherenkovs
	t->SetBranchAddress("DCPart"   ,&DCPart ); // Number of particles observed in the Drift Chambers
	t->SetBranchAddress("ECPart"   ,&ECPart ); // Number of particles observed in the ECal
	t->SetBranchAddress("SCPart"   ,&SCPart ); // Number of particles observed in the Scintillation Counters (ToFs)
	t->SetBranchAddress("Stat"     ,Stat    ); // Global status for each particle candidate
	t->SetBranchAddress("StatDC"   ,StatDC  ); // Drift Chamber status for each particle candidate
	t->SetBranchAddress("StatCC"   ,StatCC  ); // Cherenkov status for each particle
	t->SetBranchAddress("StatEC"   ,StatEC  ); // ECal status for each particle
	t->SetBranchAddress("StatSC"   ,StatSC  ); // Scintillation counter status for each particle  
	t->SetBranchAddress("particle" ,id_guess); // Guess at the particle ID made by the recon software (maybe not reliable)
	t->SetBranchAddress("EC_in"    ,EC_in   ); // Inner layer of ECal for each particle
	t->SetBranchAddress("EC_out"   ,EC_out  ); // Outer layer of ECal for each particle
	t->SetBranchAddress("EC_tot"   ,EC_tot  ); // Total energy deposit in the ECal for each particle
	t->SetBranchAddress("Nphe"     ,Nphe    ); // Number of photo-electrons per hit in the Cherenkov detectors
	t->SetBranchAddress("SC_Time"  ,SC_Time ); // Time in the scintillators per particle
	t->SetBranchAddress("SC_Path"  ,SC_Path ); // Path Length per particle
	t->SetBranchAddress("Charge"   ,Charge  ); // Charge per particle
	t->SetBranchAddress("Beta"     ,Beta    ); // Beta per particle
	t->SetBranchAddress("Mass"     ,mass    ); // Mass per particle
	t->SetBranchAddress("Momentum" ,mom     ); // Momentum magnitude per particle
	t->SetBranchAddress("Momentumx",px      ); // Momentum x component per particle
	t->SetBranchAddress("Momentumy",py      ); // Momentum y component per particle
	t->SetBranchAddress("Momentumz",pz      ); // Momentum z component per particle
	t->SetBranchAddress("Theta"    ,theta   ); // Theta per particle
	t->SetBranchAddress("Phi"      ,phi     ); // Phi per particle
	t->SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
	t->SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
	t->SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
	t->SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter

	// ---------------------------------------
	// Diagnostic electron histograms
	TH1D * h1_e_Nphe0     = new TH1D("h1_e_Nphe0"     ,"e- before cuts;# photo-electrons in CC;Counts"  ,200,   0.,200.);
	TH1D * h1_e_EC_in0    = new TH1D("h1_e_EC_in0"    ,"e- before cuts;E_in [GeV];Counts"               ,300,   0.,  1.);
	TH1D * h1_e_EC_out0   = new TH1D("h1_e_EC_out0"   ,"e- before cuts;E_out [GeV];Counts"              ,300,   0.,  1.);
	TH1D * h1_e_EC_tot0   = new TH1D("h1_e_EC_tot0"   ,"e- before cuts;E_tot [GeV];Counts"              ,300,   0.,  1.);
	TH2D * h2_e_phiTheta0 = new TH2D("h2_e_phiTheta0" ,"e- before  cuts;#phi [deg];#theta [deg];Counts" ,300,-100.,380.,300,  10., 50.);
	TH2D * h2_e_Ein_Eout0 = new TH2D("h2_e_Ein_Eout0" ,"e- before cuts;E_in/p;E_out/p;Counts"           ,300,   0., 0.5,300,   0., 0.5);
	TH2D * h2_e_Ein_Eout_0= new TH2D("h2_e_Ein_Eout_0","e- before cuts;E_in;E_out;Counts"               ,300,   0., 1.5,300,   0., 1.5);
	TH2D * h2_e_xyEC_hit0 = new TH2D("h2_e_xyEC_hit0" ,"e- before cuts;ECx [cm];ECy [cm];Counts"        ,300,-400.,400.,300,-400.,400.);
	TH2D * h2_e_p_Etot0   = new TH2D("h2_e_p_Etot0"   ,"e- before cuts;p [GeV];E_tot/p"                 ,300,   0.,  6.,300,   0., 0.7);
	TH2D * h2_e_p_E0      = new TH2D("h2_e_p_E0"      ,"e- before cuts;p [GeV];E;Counts"                ,300,   0.,  6.,300,   0.,  2.);
	TH2D * h2_e_thetaMom0 = new TH2D("e_thetaMom0"    ,"e- before cuts;#theta [deg];Mom [GeV];Counts"   ,300,  10., 50.,300,   0.,  6.);

	// Histograms to get electron PID parameters
	const float min_x = 0.5 ;
        const float max_x = 4.0 ;
	const float min   = 0.15;
	const float max   = 0.40;
	const int nSlices = 15;
	int sdiv = (int) epratio_sig_cutrange;

	TH2F * h2_e_p_Etot    = new TH2F("h2_e_p_Etot"    ,";p [GeV];E_tot/p" ,nSlices, min_x, max_x,300, min, max);
	TH1D ** h1_slices     = new TH1D * [nSlices];

	TH2D * h2_e_p_Etot1   = new TH2D("h2_e_p_Etot1"   ,Form("e- after %i#sigma pid cut;p [GeV];E_tot/p",sdiv) ,300, 0.,6.,300,0., 0.7);

	// --------------------------------------------------------------------------------------------------
	// Obtaining run number and other important parameters
	int tab_E1, tab_torus, tab_mini;
	string tab_targ;

	tab_E1    = in_num_part;
	tab_torus = 2250   ;
	tab_mini  = 5996   ;
	tab_targ  = argv[3];

	cout << "Ebeam  = " << tab_E1    << endl << "Torus  = " << tab_torus << endl;
	cout << "Mini   = " << tab_mini  << endl << "Target = " << tab_targ  << endl;

	Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, false);  // Create an instance of the Fiducial Class

	// --------------------------------------------------------------------------------------------------
        // Open up the output file
        TFile * outfile = new TFile(Form("el_Epratio_mom_%i_sim.root",tab_E1),"RECREATE");

	// --------------------------------------------------------------------------------------------------
	// Loop over events
	for (int event=0; event < nEvents ; event++)
	{
		if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}
		t->GetEvent(event);
		if (gPart <= 0) continue; // Ignore events that have no particle candidates
		double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC	

		// ---------------------------------------------------------------------------------------
		// Electron general cuts
		if (!(			(StatEC[0] > 0) && // EC status is good for the electron candidate
					(StatDC[0] > 0) && // DC status is good for the electron candidate
					(StatCC[0] > 0) && // CC status is good for the electron candidate
					(StatSC[0] > 0) && // SC status is good for the electron candidate
					(Charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
		{continue;}
		// ---------------------------------------------------------------------------------------

		h1_e_Nphe0     -> Fill( Nphe        [0]  );
		h1_e_EC_in0    -> Fill( EC_in       [0]  );
		h1_e_EC_out0   -> Fill( EC_out      [0]  );
		h1_e_EC_tot0   -> Fill( EC_tot      [0]  );
		h2_e_phiTheta0 -> Fill( phi         [0], theta        [0]);
		h2_e_Ein_Eout0 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
		h2_e_Ein_Eout_0-> Fill( EC_in[0]       , EC_out[0]       );
		h2_e_p_Etot0   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		h2_e_xyEC_hit0 -> Fill( EC_X        [0], EC_Y         [0]);
		h2_e_p_E0      -> Fill( mom         [0], el_cand_EC      );
		h2_e_thetaMom0 -> Fill( theta       [0], mom          [0]);

		h2_e_p_Etot    -> Fill( mom         [0], EC_tot[0]/mom[0]);

	}
	cerr << "Finished with the event loop...\n";

	gStyle->SetOptStat(0);

	// -------------------------------------------------------------------------------------------------- 
	// Doing fits
	TF1 * f_gaus  = new TF1("f_gaus" ,"gaus",min,max);
	TF1 * f_mean  = new TF1("f_mean" ,"pol5",min,max);
	TF1 * f_sig   = new TF1("f_sig"  ,"pol5",min,max);

	double bin_ctr [nSlices];
	double fit_mean[nSlices];
	double fit_sd  [nSlices];

	for(int i = 0 ; i < nSlices ; i++){
		bin_ctr  [i] = h2_e_p_Etot -> GetXaxis() -> GetBinCenter(i+1);
                h1_slices[i] = h2_e_p_Etot -> ProjectionY(Form("%i",i),i+1,i+2);
		h1_slices[i] -> SetTitle(Form("Bin center: p = %.2f GeV",bin_ctr[i]));

                h1_slices[i] -> Fit(f_gaus,"","",min,max);

                fit_mean[i] = f_gaus -> GetParameter(1);
                fit_sd  [i] = f_gaus -> GetParameter(2);
        }
	
	TGraphErrors * g_fit_mean = new TGraphErrors(nSlices,bin_ctr,fit_mean,0,0);
	g_fit_mean -> SetTitle("mean");
	g_fit_mean -> GetXaxis() -> SetTitle("p [GeV]");
	g_fit_mean -> GetYaxis() -> SetTitle("#mu");

	TGraphErrors * g_fit_sd   = new TGraphErrors(nSlices,bin_ctr,fit_sd  ,0,0);
	g_fit_sd -> SetTitle("standard deviation");
        g_fit_sd -> GetXaxis() -> SetTitle("p [GeV]");
        g_fit_sd -> GetYaxis() -> SetTitle("#sigma");

	// Fitting mean and standard deviation with polynomials
	g_fit_mean -> Fit(f_mean ,"","",min_x,max_x);
	g_fit_sd   -> Fit(f_sig  ,"","",min_x,max_x);
	
	// --------------------------------------------------------------------------------------------------
        // Loop over events
	double delta_t_lo_limit, delta_t_up_limit, eff_mom;

	cout << "Looping over events a second time to check the cuts" << endl;
        for (int event=0; event < nEvents ; event++)
        {
                if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}
                t->GetEvent(event);
                if (gPart <= 0) continue; // Ignore events that have no particle candidates
                double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC      

                // ---------------------------------------------------------------------------------------
                // Electron general cuts
                if (!(                  (StatEC[0] > 0) && // EC status is good for the electron candidate
                                        (StatDC[0] > 0) && // DC status is good for the electron candidate
                                        (StatCC[0] > 0) && // CC status is good for the electron candidate
                                        (StatSC[0] > 0) && // SC status is good for the electron candidate
                                        (Charge[0] < 0)    // Electron candidate curvature direction is negative
                     ))
                {continue;}
                // ---------------------------------------------------------------------------------------
		eff_mom = mom[0];
		if(eff_mom > max_x) eff_mom = max_x;

		delta_t_lo_limit = f_mean->Eval(eff_mom) - epratio_sig_cutrange * f_sig->Eval(eff_mom);
                delta_t_up_limit = f_mean->Eval(eff_mom) + epratio_sig_cutrange * f_sig->Eval(eff_mom);

		if((EC_tot[0]/mom[0]>delta_t_lo_limit)&&(EC_tot[0]/mom[0]<delta_t_up_limit)){
                	h2_e_p_Etot1   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		}
        }
        cerr << "Finished with the event loop...\n";

	// -------------------------------------------------------------------------------------------------- 
	TCanvas *c1  = new TCanvas("c1" );	h1_e_Nphe0      -> Draw();
	TCanvas *c2  = new TCanvas("c2" );	h1_e_EC_in0     -> Draw();
	TCanvas *c3  = new TCanvas("c3" );	h1_e_EC_out0    -> Draw();
	TCanvas *c4  = new TCanvas("c4" );	h1_e_EC_tot0    -> Draw();
	TCanvas *c5  = new TCanvas("c5" );	h2_e_thetaMom0  -> Draw("COLZ"); 
	TCanvas *c6  = new TCanvas("c6" );	h2_e_Ein_Eout0  -> Draw("COLZ");
	TCanvas *c7  = new TCanvas("c7" );	h2_e_Ein_Eout_0 -> Draw("COLZ");
	TCanvas *c8  = new TCanvas("c8" );	h2_e_xyEC_hit0  -> Draw("COLZ");
	TCanvas *c9  = new TCanvas("c9" );	h2_e_p_Etot0    -> Draw("COLZ");
	TCanvas *c10 = new TCanvas("c10");	h2_e_p_E0       -> Draw("COLZ");
	TCanvas *c11 = new TCanvas("c11");	h2_e_phiTheta0  -> Draw("COLZ");
	TCanvas *c12 = new TCanvas("c12");	h2_e_p_Etot     -> Draw("COLZ");

	TCanvas *c13 = new TCanvas("c13");
	c13 -> Divide(nSlices/3,3);
	for(int i = 0 ; i < nSlices ; i++){	c13 -> cd(i+1);	h1_slices[i] -> Draw();}

	TCanvas *c14 = new TCanvas("c14");
        c14 -> Divide(2,1);
	c14 -> cd(1);	g_fit_mean -> Draw("AL");
	c14 -> cd(2);	g_fit_sd   -> Draw("AL");

	TCanvas *c15 = new TCanvas("c15");
        c15 -> Divide(2,1);
        c15 -> cd(1);	h2_e_p_Etot0 -> Draw("COLZ");
	c15 -> cd(2);	h2_e_p_Etot1 -> Draw("COLZ");

	// --------------------------------------------------------------------------------------------------
	// Print histograms on a pdf file
	c1  -> Print(Form("./plots_%d.pdf(",tab_E1) ,"pdf");
	c2  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c3  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c4  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c5  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");   
	c6  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c7  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c8  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c9  -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c10 -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c11 -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c12 -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c13 -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c14 -> Print(Form("./plots_%d.pdf" ,tab_E1) ,"pdf");
	c15 -> Print(Form("./plots_%d.pdf)",tab_E1) ,"pdf");
	
	// --------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();	
	f_mean ->Write();
	f_sig  ->Write();

	// Clean up
	f->Close();
	outfile->Close();

	return 0;
}
