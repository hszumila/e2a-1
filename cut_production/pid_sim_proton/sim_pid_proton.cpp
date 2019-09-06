// ======================================================
// CODE TO DETERMINE PARAMETERS NEEDED FOR PROTON PID
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
#include <cstdlib>

#include <TGraphErrors.h>
#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\tsim_pid_proton /path/to/input/file Ebeam_MeV nuclear_target\n\n";
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
	      pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], 
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
	t->SetBranchAddress("TargetZ"  ,targetZ ); // Target Z per particle
	t->SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
	t->SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
	t->SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
	t->SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter

	// ---------------------------------------
	// Diagnostic electron histograms

        TH1D * h1_e_Nphe0    = new TH1D("h1_e_Nphe0"    ,"e- before cuts;# photo-electrons in CC;Counts"     ,200,   0.,200.);
        TH1D * h1_e_Nphe2    = new TH1D("h1_e_Nphe2"    ,"e- passing PID+fid;# photo-electrons in CC;Counts" ,200,   0.,200.);

        TH1D * h1_e_EC_in0   = new TH1D("h1_e_EC_in0"   ,"e- before cuts;E_in [GeV];Counts"                  ,300,   0.,  1.);
       	TH1D * h1_e_EC_in2   = new TH1D("h1_e_EC_in2"   ,"e- passing PID+fid;E_in [GeV];Counts"              ,300,   0.,  1.);

        TH1D * h1_e_EC_out0  = new TH1D("h1_e_EC_out0"  ,"e- before cuts;E_out [GeV];Counts"                 ,300,   0.,  1.);
        TH1D * h1_e_EC_out2  = new TH1D("h1_e_EC_out2"  ,"e- passing PID+fid;E_out [GeV];Counts"             ,300,   0.,  1.);

        TH1D * h1_e_EC_tot0   = new TH1D("h1_e_EC_tot0"  ,"e- before cuts;E_tot [GeV];Counts"                ,300,   0.,  1.);
        TH1D * h1_e_EC_tot2   = new TH1D("h1_e_EC_tot2"  ,"e- passing PID+fid;E_tot [GeV];Counts"            ,300,   0.,  1.);

        TH2D * h2_e_phiTheta0= new TH2D("h2_e_phiTheta0","e- before  cuts;#phi [deg];#theta [deg];Counts"      ,300,-100.,380.,300,10.,50.);
        TH2D * h2_e_phiTheta2= new TH2D("h2_e_phiTheta2","e- passing PID+fid;#phi [deg];#theta [deg];Counts"   ,300,-100.,380.,300,10.,50.);

        TH2D * h2_e_Ein_Eout0= new TH2D("h2_e_Ein_Eout0","e- before cuts;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
        TH2D * h2_e_Ein_Eout2= new TH2D("h2_e_Ein_Eout2","e- passing PID+fid;E_in/p;E_out/p;Counts"          ,300,   0., 0.5,300, 0.,0.5);

        TH2D * h2_e_Ein_Eout_0=new TH2D("h2_e_Ein_Eout_0","e- before cuts;E_in;E_out;Counts"                 ,300,   0., 1.5,300, 0.,1.5);
        TH2D * h2_e_Ein_Eout_2=new TH2D("h2_e_Ein_Eout_2","e- passing PID+fid;E_in;E_out;Counts"             ,300,   0., 1.5,300, 0.,1.5);

        TH2D * h2_e_xyEC_hit0= new TH2D("h2_e_xyEC_hit0","e- before cuts;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
        TH2D * h2_e_xyEC_hit2= new TH2D("h2_e_xyEC_hit2","e- passing PID+fid;ECx [cm];ECy [cm];Counts"       ,300,-400.,400.,300,-400.,400.);

        TH2D * h2_e_p_Etot0  = new TH2D("h2_e_p_Etot0"  ,"e- before cuts;p [GeV];E_tot/p;Counts"             ,300,   0.,  7.,300, 0.,0.7);
        TH2D * h2_e_p_Etot2  = new TH2D("h2_e_p_Etot2"  ,"e- passing PID+fid;p [GeV];E_tot/p;Counts"         ,300,   0.,  7.,300, 0.,0.7);

        TH2D * h2_e_p_E0     = new TH2D("h2_e_p_E0"     ,"e- before cuts;p [GeV];E;Counts"                   ,300,   0.,  7.,300, 0.,2.);
        TH2D * h2_e_p_E2     = new TH2D("h2_e_p_E2"     ,"e- before PIF+fid;p [GeV];E;Counts"                ,300,   0.,  7.,300, 0.,2.);

        TH2D * h2_e_thetaMom0= new TH2D("e_thetaMom0"     ,"e- before cuts;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
        TH2D * h2_e_thetaMom2= new TH2D("e_thetaMom2"     ,"e- passing PID+fid;#theta [deg];Mom [GeV];Counts"   ,300,  10., 50.,300, 0., 6.);

        // ---------------------------------------
        // Diagnostic positive particle histograms
        TH1D * h1_p_mass     = new TH1D("h1_pos_mass"   ,"+  passing fid. cuts;mass [GeV];Counts"          ,300,  0.,3.5 );
        TH2D * h2_p_pMass    = new TH2D("h2_pos_pMass"  ,"+  passing fid. cuts;p [GeV];mass [GeV];Counts"  ,300,  0., 4.,300, 0.,3.5);
        TH2D * h2_pos_pBeta  = new TH2D("h2_pos_pBeta"  ,"+  passing fid. cuts;p [GeV];#beta;Counts"       ,300,  0., 4.,300, 0.,1.3);

        // ---------------------------------------
        // Diagnostic proton histograms
	TH2D * h2_p_phiTheta0= new TH2D("h2_p_phiTheta0","p before cuts;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
        TH2D * h2_p_deltaTmom0=new TH2D("h2_p_deltaTmom0","p before cuts;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 6.);

        TH2D * h2_p_pBeta    = new TH2D("h2_p_pBeta"    ,"p passing fid. cuts;p [GeV];#beta;Counts"                 ,300,   0.,  4.,300, 0.,1.3);

	// ---------------------------------------
	// Histograms to get proton PID parameters
	const float min_x = 0.5;
        const float max_x = 3.3;
	const float min   = -1.;
	const float max   =  1.;
	const int nSlices =  12;
	int sdiv = (int) pdeltat_sig_cutrange;

	TH2F * h2_p_mom_deltaT0    = new TH2F("h2_p_mom_deltaT0"    ,";p [GeV];#Delta t [ns]" ,nSlices, min_x, max_x,50, min, max);
	TH1D ** h1_slices          = new TH1D * [nSlices];

	TH2D * h2_p_mom_deltaT1= new TH2D("h2_p_mom_deltaT1",     "p+ before pid cut;p [GeV];#Delta t [ns]"               ,300,  0.,  5.,300, -4., 4.);
	TH2D * h2_p_mom_deltaT2= new TH2D("h2_p_mom_deltaT2",Form("p+ after %i#sigma pid cut;p [GeV];#Delta t [ns]",sdiv) ,300,  0.,  5.,300, -4., 4.);
	
	TH2D * h2_prot_pBeta1  = new TH2D("h2_prot_pBeta1"  ,     "p+ before pid cut;p [GeV];#beta"                       ,300,  0., 4.,300, 0.,1.3);
	TH2D * h2_prot_pBeta2  = new TH2D("h2_prot_pBeta2"  ,Form("p+ after %i#sigma pid cut;p [GeV];#beta",sdiv)         ,300,  0., 4.,300, 0.,1.3);

	// ---------------------------------------
	// Setting up output tree and branches
	double e_vz_corrected; //, e_vz, e_mom[3], e_phi_mod;
	double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
	double EC_in_cut, el_EC_cut;
	double e_t0,beta_assuming_proton,p_t0,delta_t;

	TVector3 e_ec_xyz;
	TVector3 T3_e_mom, T3_p_mom;

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
        TFile * outfile = new TFile(Form("protdeltat_mom_%i_%i_sim.root",tab_E1,tab_torus),"RECREATE");

	// Values for some cuts
	if      (tab_E1 == 4461){
		EC_in_cut = 0.055; //GeV (Values for Energy deposited in EC inner layer cut)
		el_EC_cut = 0.330; //GeV (Values for Enough total energy in the EC cut)
	}
	else if (tab_E1 == 2261){
		EC_in_cut = 0.060; //GeV (Values for Energy deposited in EC inner layer cut)
		el_EC_cut = -999.; // No cut in this case (Values for Enough total energy in the EC cut)
	}
	else {
		cout << "Error: Check skim_tree and add parameters for Ebeam = " << tab_E1 << endl;
		exit(-2);
	}

	// --------------------------------------------------------------------------------------------------
	// Loop over events
	for (int event=0; event < nEvents ; event++)
	{
		if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}

		t->GetEvent(event);	

		if (gPart <= 0) continue; // Ignore events that have no particle candidates	

		// --------------------------------------------------------------------------------------------------
		// Sector index for electrons
		int e_sect = (int)(phi[0]+30)/60;
		if (e_sect>5) e_sect = 5;
		if (e_sect<0) e_sect = 0;

		// --------------------------------------------------------------------------------------------------
		double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC
		T3_e_mom.SetXYZ(px[0],py[0],pz[0]); // Electron momentum expressed in a TVector3
		e_vz_corrected = targetZ[0]+fid_params.vz_corr(T3_e_mom);
		e_ec_xyz.SetXYZ(EC_X[0],EC_Y[0],EC_Z[0]);

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

		// ---------------------------------------------------------------------------------------
		//Electron particle Identification
		if (!(                  (EC_in [0] > EC_in_cut) &&      // Electron candidate has enough energy deposit in inner layer of EC    
					(el_cand_EC > el_EC_cut) &&     // Enough total energy in the EC
					(fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) // Electron PID (E/p)
		     ))
		{continue;}

		// ---------------------------------------
		// Additional cut for 2GeV data:
		/*
		double el_sccc_dt = SC_Time[0] - CC_Time[0] - (SC_Path[0] - CC_Path[0])/(c_m_s*ns_to_s*100.);

		if(			(tab_E1==2261)&&(	
					CC_Chi2[0]>=0.1 ||
					el_sccc_dt < sc_cc_dt_cut_sect[e_sect] ||
					sqrt(mom[0]*mom[0]+me*me)>tab_E1/1000.
					))
		{continue;}	
		*/

		// ---------------------------------------------------------------------------------------
		// Electron Fiducial cuts
		if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
		if (!fid_params.CutUVW_e(e_ec_xyz)       ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);

		// Fill some diagnostic histograms
		h1_e_Nphe2     -> Fill( Nphe        [0]  );
                h1_e_EC_in2    -> Fill( EC_in       [0]  );
                h1_e_EC_out2   -> Fill( EC_out      [0]  );
                h1_e_EC_tot2   -> Fill( EC_tot      [0]  );
                h2_e_phiTheta2 -> Fill( phi         [0], theta        [0]);
                h2_e_Ein_Eout2 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
                h2_e_Ein_Eout_2-> Fill( EC_in[0]       , EC_out[0]       );
                h2_e_p_Etot2   -> Fill( mom         [0], EC_tot[0]/mom[0]);
                h2_e_xyEC_hit2 -> Fill( EC_X        [0], EC_Y         [0]);
                h2_e_p_E2      -> Fill( mom         [0], el_cand_EC      );
                h2_e_thetaMom2 -> Fill( theta       [0], mom          [0]);
		
		// --------------------------------------------------------------------------------------------------
		// Loop over events looking for other particles
		for (int i=1 ; i<gPart ; i++)
		{	
			T3_p_mom.SetXYZ(px[i],py[i],pz[i]);	

			e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;

			beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
			p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
			delta_t = p_t0 - e_t0;

			// ------------------------------------------------------------------------------------------
			// Test if positive particle
			if(             (StatSC[i] > 0) && 		// SC status is good for the positive candidate
					(StatDC[i] > 0) &&              // DC status is good for the positive candidate
					(Stat  [i] > 0) &&		// Global status is good for the positive candidate
					(Charge[i] > 0) 		// Charge is positive
			  )
			{

				// Passing positive hadron fiducial cuts
				if(fid_params.pFiducialCut(T3_p_mom)){

					h1_p_mass         -> Fill(mass[i]);
					h2_p_phiTheta0    -> Fill(phi [i]    ,theta[i]);
					h2_p_pMass        -> Fill(mom [i]    ,mass [i]);
					h2_pos_pBeta      -> Fill(mom [i]    ,Beta [i]);
					h2_p_deltaTmom0   -> Fill(delta_t    ,mom  [i]);
					
					h2_p_mom_deltaT0  -> Fill(mom [i]    ,delta_t);

				}
			}
		}
	}
	cerr << "Finished with the event loop...\n";

	gStyle->SetOptStat(0);

	// -------------------------------------------------------------------------------------------------- 
	// Doing fits
	TF1 * f_gaus  = new TF1("f_gaus" ,"gaus",min,max);
	TF1 * mean_pol9  = new TF1("mean_pol9" ,"pol5",min,max);
	TF1 * sig_pol9   = new TF1("sig_pol9"  ,"pol5",min,max);

	double bin_ctr [nSlices];
	double fit_mean[nSlices];
	double fit_sd  [nSlices];

	for(int i = 0 ; i < nSlices ; i++){
		bin_ctr  [i] = h2_p_mom_deltaT0 -> GetXaxis() -> GetBinCenter(i+1);
                h1_slices[i] = h2_p_mom_deltaT0 -> ProjectionY(Form("%i",i),i+1,i+2);
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
	g_fit_mean -> Fit(mean_pol9 ,"","",min_x,max_x);
	g_fit_sd   -> Fit(sig_pol9  ,"","",min_x,max_x);

	// --------------------------------------------------------------------------------------------------
        // Loop over events

	double prot_mom_lim = 2.;
        double prot_min = 0.3; //GeV

	double delta_t_up_limit, delta_t_lo_limit, mom_eff;

	cout << "Looping over events a second time to check the cuts" << endl;
	for (int event=0; event < nEvents ; event++)
        {
                if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}

                t->GetEvent(event);

                if (gPart <= 0) continue; // Ignore events that have no particle candidates     

                // --------------------------------------------------------------------------------------------------
                // Sector index for electrons
                int e_sect = (int)(phi[0]+30)/60;
                if (e_sect>5) e_sect = 5;
                if (e_sect<0) e_sect = 0;

                // --------------------------------------------------------------------------------------------------
                double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC
                T3_e_mom.SetXYZ(px[0],py[0],pz[0]); // Electron momentum expressed in a TVector3
                e_vz_corrected = targetZ[0]+fid_params.vz_corr(T3_e_mom);
                e_ec_xyz.SetXYZ(EC_X[0],EC_Y[0],EC_Z[0]);

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
                //Electron particle Identification
		if (!(                  (EC_in [0] > EC_in_cut) &&      // Electron candidate has enough energy deposit in inner layer of EC    
                                        (el_cand_EC > el_EC_cut) &&     // Enough total energy in the EC
                                        (fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) // Electron PID (E/p)
                     ))
                {continue;}

                // ---------------------------------------
                // Additional cut for 2GeV data:
                /*
                double el_sccc_dt = SC_Time[0] - CC_Time[0] - (SC_Path[0] - CC_Path[0])/(c_m_s*ns_to_s*100.);

                if(                     (tab_E1==2261)&&(       
                                        CC_Chi2[0]>=0.1 ||
                                        el_sccc_dt < sc_cc_dt_cut_sect[e_sect] ||
                                        sqrt(mom[0]*mom[0]+me*me)>tab_E1/1000.
                                        ))
                {continue;}     
                */

                // ---------------------------------------------------------------------------------------
                // Electron Fiducial cuts
                if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
                if (!fid_params.CutUVW_e(e_ec_xyz)     ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);

                // --------------------------------------------------------------------------------------------------
                // Loop over events looking for other particles
                for (int i=1 ; i<gPart ; i++)
                {
                        T3_p_mom.SetXYZ(px[i],py[i],pz[i]); 

                        e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;

                        beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
                        p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
                        delta_t = p_t0 - e_t0;

                        // ------------------------------------------------------------------------------------------
                        // Test if positive particle
			if(             (StatSC[i] > 0) &&              // SC status is good for the positive candidate
                                        (StatDC[i] > 0) &&              // DC status is good for the positive candidate
                                        (Stat  [i] > 0) &&              // Global status is good for the positive candidate
                                        (Charge[i] > 0)                 // Charge is positive
                          )
                        {

                                // Passing positive hadron fiducial cuts
                                if(fid_params.pFiducialCut(T3_p_mom)){

					h2_p_mom_deltaT1  -> Fill(mom [i]    ,delta_t );
					h2_prot_pBeta1    -> Fill(mom [i]    ,Beta [i]);

					mom_eff = mom[i];
					if (mom_eff > prot_mom_lim) mom_eff = prot_mom_lim;
        				
					if (mom_eff > prot_min    ){ 
        					delta_t_up_limit = mean_pol9->Eval(mom_eff) + pdeltat_sig_cutrange * sig_pol9->Eval(mom_eff);
        					delta_t_lo_limit = mean_pol9->Eval(mom_eff) - pdeltat_sig_cutrange * sig_pol9->Eval(mom_eff);

						if((delta_t<delta_t_up_limit)&&(delta_t>delta_t_lo_limit)){
							h2_p_mom_deltaT2  -> Fill(mom [i]    ,delta_t );
							h2_prot_pBeta2    -> Fill(mom [i]    ,Beta [i]);
						}
					}
                                }
                        }
                }
        }
        cerr << "Finished with the event loop...\n";

	// -------------------------------------------------------------------------------------------------- 
	TCanvas *c1  = new TCanvas("c1" );
	c1 -> Divide(2,2);
	c1 -> cd(1);	h1_e_Nphe0      -> Draw();	h1_e_Nphe2      -> Draw("same");
	c1 -> cd(2);	h1_e_EC_in0     -> Draw();	h1_e_EC_in2     -> Draw("same");
	c1 -> cd(3);	h1_e_EC_out0    -> Draw();	h1_e_EC_out2    -> Draw("same");
	c1 -> cd(4);	h1_e_EC_tot0    -> Draw();	h1_e_EC_tot2    -> Draw("same");

	TCanvas *c2  = new TCanvas("c2" );
        c2 -> Divide(2,2);
	c2 -> cd(1);	h2_e_Ein_Eout0  -> Draw("COLZ");
	c2 -> cd(2);    h2_e_Ein_Eout2  -> Draw("COLZ");
	c2 -> cd(3);	h2_e_Ein_Eout_0 -> Draw("COLZ");
	c2 -> cd(4);    h2_e_Ein_Eout_2 -> Draw("COLZ");

	TCanvas *c3  = new TCanvas("c3" );
        c3 -> Divide(2,2);
	c3 -> cd(1);	h2_e_thetaMom0  -> Draw("COLZ");
	c3 -> cd(2);    h2_e_thetaMom2  -> Draw("COLZ"); 
	c3 -> cd(3);	h2_e_xyEC_hit0  -> Draw("COLZ");
	c3 -> cd(4);    h2_e_xyEC_hit2  -> Draw("COLZ");
	
	TCanvas *c4  = new TCanvas("c4" );
        c4 -> Divide(2,2);
	c4 -> cd(1);	h2_e_p_Etot0    -> Draw("COLZ");
	c4 -> cd(2);    h2_e_p_Etot2    -> Draw("COLZ");
	c4 -> cd(3);	h2_e_p_E0       -> Draw("COLZ");
	c4 -> cd(4);    h2_e_p_E2       -> Draw("COLZ");

	TCanvas *c5  = new TCanvas("c5" );
        c5 -> Divide(1,2);
	c5 -> cd(1);	h2_e_phiTheta0  -> Draw("COLZ");
	c5 -> cd(2);    h2_e_phiTheta2  -> Draw("COLZ");

	TCanvas *c6   = new TCanvas("c6"  );	h1_p_mass       -> Draw();
        TCanvas *c7   = new TCanvas("c7"  );	h2_p_phiTheta0  -> Draw("COLZ");
        TCanvas *c8   = new TCanvas("c8"  );	h2_p_pMass      -> Draw("COLZ");
        TCanvas *c9   = new TCanvas("c9"  );	h2_pos_pBeta    -> Draw("COLZ");
        TCanvas *c10  = new TCanvas("c10" );	h2_p_deltaTmom0 -> Draw("COLZ");

	TCanvas *c11  = new TCanvas("c11" );    h2_p_mom_deltaT0-> Draw("COLZ");

	TCanvas *c12 = new TCanvas("c12");
	c12 -> Divide(nSlices/3,3);
	for(int i = 0 ; i < nSlices ; i++){	c12 -> cd(i+1);	h1_slices[i] -> Draw();}

	TCanvas *c13 = new TCanvas("c13");
        c13 -> Divide(2,1);
	c13 -> cd(1);	g_fit_mean -> Draw("AL");
	c13 -> cd(2);	g_fit_sd   -> Draw("AL");

	TCanvas *c14 = new TCanvas("c14");
        c14 -> Divide(2,1);
        c14 -> cd(1);	h2_p_mom_deltaT1 -> Draw("COLZ");
	c14 -> cd(2);	h2_p_mom_deltaT2 -> Draw("COLZ");
	
	TCanvas *c15 = new TCanvas("c15");
        c15 -> Divide(2,1);
	c15 -> cd(1);	h2_prot_pBeta1 -> Draw("COLZ");
	c15 -> cd(2);	h2_prot_pBeta2 -> Draw("COLZ");

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
	// Write the output file
	outfile->cd();	
	mean_pol9 ->Write();
	sig_pol9  ->Write();

	// Clean up
	f->Close();
	outfile->Close();
	
	return 0;
}
