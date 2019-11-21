#include "Fiducial.h"
#include "e2a_constants.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TFile.h"
// ===================================================================================================================================
Fiducial::Fiducial(int E_beam, int torus, int mini, std::string target, bool data)
{
  // Initialize some variables
  el_Ep_ratio_mean = NULL;
  el_Ep_ratio_sig  = NULL;
  prot_deltat_sig  = NULL;
  prot_deltat_mean = NULL;
  pip_deltat_sig   = NULL;
  pip_deltat_mean  = NULL;
  pim_deltat_sig   = NULL;
  pim_deltat_mean  = NULL;
  
  // Initialize the key run settings
  E1=E_beam;
  torus_current = torus;
  mini_current = mini;
  tar = target;
  is_data = data;
  
  e2adir = std::string(getenv("E2A_INSTALL"));
    
  // Read in the various parameters
  bool all_ok=true;
  all_ok &= read_e_fid_params      (); // Electron fiducial regions
  std::cerr << "after e_fid" << " " << all_ok << "\n";
  all_ok &= read_e_pcor_params     (); // Electron momentum corrections
  std::cerr << "after e_pcor" << " " << all_ok << "\n";
  all_ok &= read_e_pid_params      (); // Electron E/p params
  std::cerr << "after e_pid" << " " << all_ok << "\n";
  all_ok &= read_p_pid_params      (); // Proton delta_t vs mom params
  std::cerr << "after p_pid" << " " << all_ok << "\n";
  all_ok &= read_pip_pid_params    (); // Pi+ delta_t vs mom params
  std::cerr << "after pip_pid" << " " << all_ok << "\n";
  all_ok &= read_pim_pid_params    (); // Pi- delta_t vs mom params
  std::cerr << "after pim_pid" << " " << all_ok << "\n";
  all_ok &= read_vz_cor_params     (); // vz corrections
  std::cerr << "after vz_cor" << " " << all_ok << "\n";
  all_ok &= read_p_fid_params      (); // Proton fiducial regions
  std::cerr << "after p_fid" << " " << all_ok << "\n";
  all_ok &= read_pim_fid_params      (); // Piminus fiducial regions
  std::cerr << "after pim_fid" << " " << all_ok << "\n";
  all_ok &= read_n_pathlength_corr (); // Neutron pathlength correction params
  std::cerr << "after n_path" << " " << all_ok << "\n";
  if (all_ok)
    std::cerr << "Successfully read in the various parameters...\n";
  else
    {
      std::cerr << "Failed to read in the parameter files. Exiting...\n";
      exit(-1);
    }  
}
// ===================================================================================================================================
Fiducial::~Fiducial()
{
	// Memory clean up
	if (prot_deltat_sig )	delete prot_deltat_sig;
	if (prot_deltat_mean)	delete prot_deltat_mean;
	if (pip_deltat_sig  )   delete pip_deltat_sig;
  if (pip_deltat_mean )   delete pip_deltat_mean;
	if (pim_deltat_sig  )   delete pim_deltat_sig;
  if (pim_deltat_mean )   delete pim_deltat_mean;
	if (el_Ep_ratio_sig )	delete el_Ep_ratio_sig;
	if (el_Ep_ratio_mean)	delete el_Ep_ratio_mean;
}
// ===================================================================================================================================
bool Fiducial::pim_inFidRegion(TVector3 mom)
{
  // Assemble sanitized phi and sector variables
  double phi = mom.Phi();
  if (phi < -M_PI/6.) phi+= 2.*M_PI;
  int sector = (phi+M_PI/6.)/(M_PI/3.);
  sector = sector%6;
  double phi_deg = phi * 180./M_PI; 
  phi_deg -= sector*60; // phi value within sector

  double theta = mom.Theta();
  double theta_deg = theta * 180./M_PI;
  double mom_pi = mom.Mag();

  // Memory for assembling the parameters
  Double_t phipars[5]={0,0,0,0,0};

  // We need to add santization bounds for the different beam/torus combos, but for now...
  if (mom_pi < .15)
    mom_pi = .15;
  if (mom_pi > 1.1)
    mom_pi = 1.1;

  // Assemble the parameters
  for(Int_t mompar=0;mompar<6;mompar++)
    for(Int_t phipar=0;phipar<5;phipar++)
        phipars[phipar]+=fgPar_Pimfid[sector][phipar][mompar]*pow(mom_pi,mompar);

  // Test if our theta is more forward than the "pinch"
  if (theta_deg < phipars[4])
    return false;

  // Test if we are within the phi bounds
  Int_t uplow;
  Double_t phicutoff;
  if(phi_deg<=0) // If in the lower-phi half of the sector
    {
      phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
      if (phi_deg < phicutoff)
	return false;
    }
  else
    {
      phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
      if (phi_deg > phicutoff)
	return false;
    }

  // Now we have to evaluate theta gaps. It looks like these are only for 1GeV running so far...
  bool SCpdcut = true; // If we actually want to evaluate gaps
  bool status = true;

  if (SCpdcut) // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
    {
      if (sector == 2){
	float badpar3[4]; // 4 parameters to determine the positions of the two theta gaps
	for (int i=0; i<4; i++){
	  badpar3[i] = 0;
	  // calculate the parameters using pol7
	  for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_pi + fgPar_Pimfid_Theta_S3[i][d];}
	}
	for(int ipar=0;ipar<2;ipar++)
	  status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
      }
      // sector 3 has one bad paddle
      else if (sector == 3){
	float badpar4[2];     // 2 parameters to determine the position of the theta gap
	for (int i=0; i<2; i++){
	  badpar4[i] = 0;
	  // calculate the parameters using pol7
	  for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_pi + fgPar_Pimfid_Theta_S4[i][d];}
	}
	status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
      }
      // sector 4 has four bad paddles
      else if (sector == 4){ 
	Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar5[i] = 0;
	  // calculate the parameters using pol7
	  for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_pi + fgPar_Pimfid_Theta_S5[i][d];}
	}
	if (mom_pi<1.25) badpar5[0] = 23.4*torus_current/2250.;
	if (mom_pi<1.27) badpar5[1] = 24.0*torus_current/2250.; // some dummy constants. see fiducial cuts webpage.
	for(Int_t ipar=0;ipar<4;ipar++)
	  status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
      }
    }
  return status;
}

bool Fiducial::e_inFidRegion(TVector3 mom)
{
  // Establish the sector;
  double phi = mom.Phi();
  if (phi < -M_PI/6.) phi+= 2.*M_PI;
  int sector = (phi+M_PI/6.)/(M_PI/3.);
  double phi_deg = phi * 180./M_PI;

  double theta = mom.Theta();
  double theta_deg = theta * 180./M_PI;
  double mom_e = mom.Mag();

  // ---------------------------------------------------------------------------
  // Cut for Ebeam = 4.4GeV and 2250A data.
  if ( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.){
    bool status = true;
    double phiMin, phiMax;
    // Sanitize theta
    if (theta_deg < 15.) return false;
    
    // Sanitize momentum
    if (mom_e > 3.7) mom_e = 3.7;
    if (mom_e < 0.9)
      {
	return false;
      }
    
    // Assemble the polynomials
    double t0 = fgPar_Efid_t0_p[sector][0]/pow(mom_e, fgPar_Efid_t0_p[sector][1]);
    double t1 = 0.; 
    double b[2]={0.,0.};
    double a[2]={0.,0.};
    for(int k=0; k<6; k++)
      {
	double mom_to_the_k = pow(mom_e,k);
	t1 += fgPar_Efid_t1_p[sector][k]*mom_to_the_k;
	for(int l=0; l<2; l++)
	  {
	    a[l] += fgPar_Efid_a_p[sector][l][k]*mom_to_the_k;
	    b[l] += fgPar_Efid_b_p[sector][l][k]*mom_to_the_k;
	  }
      }
    
    // Calculate the limits
    phiMin=sector*M_PI/3.; // Default is the center line of each sector
    phiMax=sector*M_PI/3.;
    if(t1 < 45.) t1 = 45.;
    if((t0 < theta_deg) && (theta_deg < t1))
      {
	phiMin -= M_PI/180.*b[0]*(1. - 1/((theta_deg - t0)/(b[0]/a[0]) + 1.));
	phiMax += M_PI/180.*b[1]*(1. - 1/((theta_deg - t0)/(b[1]/a[1]) + 1.));
      }
    if(mom_e <= 2.0)
      {
	bool SCpdcut = true;
	if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	  if (status){
	    int tsector = sector + 1;
	    // sector 3 has two bad paddles
	    if (tsector == 3){
	      float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	      for (int i=0; i<4; i++){
		badpar3[i] = 0;
		// calculate the parameters using pol7
		for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_Efid_Theta_S3[i][d];}
	      }
	      for(int ipar=0;ipar<2;ipar++)
		status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
	    }
	    // sector 4 has one bad paddle
	    else if (tsector == 4){
	      float badpar4[2];     // 2 parameters to determine the position of the theta gap
	      for (int i=0; i<2; i++){
		badpar4[i] = 0;
		// calculate the parameters using pol7
		for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_Efid_Theta_S4[i][d];}
	      }
	      status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
	    }
	    // sector 5 has four bad paddles
	    else if (tsector == 5){ 
	      Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	      for (Int_t i=0; i<8; i++){
		badpar5[i] = 0;
		// calculate the parameters using pol7
		for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_Efid_Theta_S5[i][d];}
	      }
	      if (mom_e<1.25) badpar5[0] = 23.4;
	      if (mom_e<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	      for(Int_t ipar=0;ipar<4;ipar++)
		status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
	    }
	  }
	}
	return (status && (phi < phiMax) && (phi>phiMin));
      }
    else{
        bool SCpdcut = true;
        if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
          if (status){
            int tsector = sector + 1;
            // sector 3 has two bad paddles
            if (tsector == 3){
              float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
              for (int i=0; i<4; i++){
                badpar3[i] = 0;
                // calculate the parameters using 1/p
                badpar3[i] = fgPar_Efid_Theta_S3_extra[i][0] + fgPar_Efid_Theta_S3_extra[i][1]/mom_e + fgPar_Efid_Theta_S3_extra[i][2]/(mom_e*mom_e) + fgPar_Efid_Theta_S3_extra[i][3]/(mom_e*mom_e*mom_e);
              }
              for(int ipar=0;ipar<2;ipar++)
                status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
            }
            // sector 4 has one bad paddle
            else if (tsector == 4){
              float badpar4[2];     // 2 parameters to determine the position of the theta gap
              for (int i=0; i<2; i++){
                badpar4[i] = 0;
                // calculate the parameters using 1/p
                badpar4[i] = fgPar_Efid_Theta_S4_extra[i][0] + fgPar_Efid_Theta_S4_extra[i][1]/mom_e + fgPar_Efid_Theta_S4_extra[i][2]/(mom_e*mom_e) + fgPar_Efid_Theta_S4_extra[i][3]/(mom_e*mom_e*mom_e);
              }
              status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
            }
            // sector 5 has four bad paddles
            else if (tsector == 5){ 
              Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
              for (Int_t i=0; i<8; i++){
                badpar5[i] = 0;
                // calculate the parameters using 1/p
                badpar5[i] = fgPar_Efid_Theta_S5_extra[i][0] + fgPar_Efid_Theta_S5_extra[i][1]/mom_e + fgPar_Efid_Theta_S5_extra[i][2]/(mom_e*mom_e) + fgPar_Efid_Theta_S5_extra[i][3]/(mom_e*mom_e*mom_e);
              }
              if (mom_e<1.25) badpar5[0] = 23.4;
              if (mom_e<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
              for(Int_t ipar=0;ipar<4;ipar++)
                status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
            }
          }
        }
        return (status && (phi < phiMax) && (phi>phiMin));
      }
  }
	// ---------------------------------------------------------------------------
	// Cut for Ebeam = 2.2GeV and 2250A data.

	else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.){

	  // Sanitize momentum
	  if (mom_e > 2.0) mom_e = 2.0;
    if (mom_e < 0.3) return false;
		bool status = true;
		phi_deg -= sector*60;
		Float_t par[6];               // six parameters to determine the outline of Theta vs Phi
		for (int i=0; i<6; i++){
			par[i] = 0;
			// calculate the parameters using pol8
			for (int d=8; d>=0; d--){par[i] = par[i]*mom_e + fgPar_Efid[sector][i][d];}
		}
		if (phi_deg < 0) {
			float tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi_deg);
			status = (theta_deg>tmptheta && tmptheta>=par[0] && theta_deg<par[1]);
		}
		else {
			float tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi_deg);
			status = (theta_deg>tmptheta && tmptheta>=par[0] && theta_deg<par[1]);
		}
		// --------
		// by now, we have checked if the electron is within the outline of theta vs phi plot
		bool SCpdcut = true;
		if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
			if (status){
				int tsector = sector + 1;
				// sector 3 has two bad paddles
				if (tsector == 3){
					float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
					for (int i=0; i<4; i++){
						badpar3[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_Efid_Theta_S3[i][d];}
					}
					for(int ipar=0;ipar<2;ipar++)
						status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
				}
				// sector 4 has one bad paddle
				else if (tsector == 4){
					float badpar4[2];     // 2 parameters to determine the position of the theta gap
					for (int i=0; i<2; i++){
						badpar4[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_Efid_Theta_S4[i][d];}
					}
					status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
				}
				// sector 5 has four bad paddles
				else if (tsector == 5){ 
					Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
					for (Int_t i=0; i<8; i++){
						badpar5[i] = 0;
						// calculate the parameters using pol7
						for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_Efid_Theta_S5[i][d];}
					}
					if (mom_e<1.25) badpar5[0] = 23.4;
					if (mom_e<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
					for(Int_t ipar=0;ipar<4;ipar++)
						status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  else if ( E1 > 1000 && E1 < 2000 && torus_current < 1510 && torus_current > 1490){
    phi_deg-=sector*60.;
    Double_t elmom = (mom.Mag())*1000;
    Double_t thetapars[5]={0,0,0,0,0};
    bool status = true;
    if (elmom < 300)
      return false;
    if (elmom > 1100)
      elmom = 1100;

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t thetapar=0;thetapar<5;thetapar++) {
        thetapars[thetapar]+=fgPar_1gev_1500_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
      }
    }
    Int_t uplow;
    Double_t thetacutoff;
    if(phi_deg<=0) {
      uplow=1;
      thetacutoff=((phi_deg*(thetapars[0]-(thetapars[1]/thetapars[2])))+
                   (double(uplow)*thetapars[2]*thetapars[0]))/(phi_deg+(double(uplow)*thetapars[2]));
    }
    else {
      uplow=-1;
      thetacutoff=( (phi_deg*(thetapars[0]-(thetapars[3]/thetapars[4]))) +
                    (double(uplow)*thetapars[4]*thetapars[0]))/(phi_deg+(double(uplow)*thetapars[4]) );
    }
    status = ((theta_deg>thetacutoff) && (thetacutoff>=thetapars[0]));
    // by now, we have checked if the electron is within the outline of theta vs phi plot
		bool SCpdcut = true;
		if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
			if (status){
				int tsector = sector + 1;
				// sector 3 has two bad paddles
				if (tsector == 3){
					float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
					for (int i=0; i<4; i++){
						badpar3[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_1gev_1500_Efid_Theta_S3[i][d];}
					}
					for(int ipar=0;ipar<2;ipar++)
						status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
				}
				// sector 4 has one bad paddle
				else if (tsector == 4){
					float badpar4[2];     // 2 parameters to determine the position of the theta gap
					for (int i=0; i<2; i++){
						badpar4[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_1gev_1500_Efid_Theta_S4[i][d];}
					}
					status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
				}
				// sector 5 has four bad paddles
				else if (tsector == 5){ 
					Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
					for (Int_t i=0; i<8; i++){
						badpar5[i] = 0;
						// calculate the parameters using pol7
						for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_1gev_1500_Efid_Theta_S5[i][d];}
					}
					if (mom_e<1.25) badpar5[0] = 23.4*1500/2250;
					if (mom_e<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
					for(Int_t ipar=0;ipar<4;ipar++)
						status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  else if ( E1 > 1000 && E1 < 2000 && torus_current < 760 && torus_current > 740){
    phi_deg-=sector*60.;
    Double_t elmom = (mom.Mag())*1000;
    Double_t thetapars[5]={0,0,0,0,0};
    bool status = true;
    if (elmom < 300)
      return false;
    if (elmom > 1100)
      elmom = 1100;

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t thetapar=0;thetapar<5;thetapar++) {
        thetapars[thetapar]+=fgPar_1gev_750_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
      }
    }
    Int_t uplow;
    Double_t thetacutoff;
    if(phi_deg<=0) {
      uplow=1;
      thetacutoff=((phi_deg*(thetapars[0]-(thetapars[1]/thetapars[2])))+
                   (double(uplow)*thetapars[2]*thetapars[0]))/(phi_deg+(double(uplow)*thetapars[2]));
    }
    else {
      uplow=-1;
      thetacutoff=( (phi_deg*(thetapars[0]-(thetapars[3]/thetapars[4]))) +
                    (double(uplow)*thetapars[4]*thetapars[0]))/(phi_deg+(double(uplow)*thetapars[4]) );
    }
    status = ((theta_deg>thetacutoff) && (thetacutoff>=thetapars[0]));
// by now, we have checked if the electron is within the outline of theta vs phi plot
		bool SCpdcut = true;
		if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
			if (status){
				int tsector = sector + 1;
				// sector 3 has two bad paddles
				if (tsector == 3){
					float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
					for (int i=0; i<4; i++){
						badpar3[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_1gev_750_Efid_Theta_S3[i][d];}
					}
					for(int ipar=0;ipar<2;ipar++)
						status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
				}
				// sector 4 has one bad paddle
				else if (tsector == 4){
					float badpar4[2];     // 2 parameters to determine the position of the theta gap
					for (int i=0; i<2; i++){
						badpar4[i] = 0;
						// calculate the parameters using pol7
						for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_1gev_750_Efid_Theta_S4[i][d];}
					}
					status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
				}
				// sector 5 has four bad paddles
				else if (tsector == 5){ 
					Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
					for (Int_t i=0; i<8; i++){
						badpar5[i] = 0;
						// calculate the parameters using pol7
						for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_1gev_750_Efid_Theta_S5[i][d];}
					}
					if (mom_e<1.25) badpar5[0] = 23.4*750/2250;
					if (mom_e<1.27) badpar5[1] = 24.0*750/2250; // some dummy constants. see fiducial cuts webpage.
					for(Int_t ipar=0;ipar<4;ipar++)
						status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  // ---------------------------------------------------------------------------
  else {
    std::cerr << "e_inFidRegion doesn't have correction parameters for the given input. Check it and fix it!\n";
    exit(-3);
  }
}
// ===================================================================================================================================
bool Fiducial::read_n_pathlength_corr()
{
  //Parameters for neutron path length correction

  char param_file_name[256];
  sprintf(param_file_name,"%s/n_pathlength_corr_%d.dat",e2adir.c_str(),E1);
  std::ifstream param_file(param_file_name);
  std::cerr<<param_file_name<<std::endl;
       
	param_file >> pl_corr_in  ;
	param_file >> pl_corr_out ;
	param_file >> pl_corr_both;
       
	param_file.close();

	if(pl_corr_in==0&&pl_corr_out==0&&pl_corr_both==0)
		std::cerr << "*** WARNING *** Won't be correcting neutron path length since there are no available parameters!" << std::endl;

  return true;
}
// ===================================================================================================================================
bool Fiducial::read_p_fid_params()
{
	//Parameters for 4 GeV proton's Fiducial Cut Rustam Niyazov
	//"http://www.physics.odu.edu/~rust/clas/fidp.html"

	char param_file_name[256];
	sprintf(param_file_name,"%s/PFID_%d_%d.dat",e2adir.c_str(),E1,torus_current);
	std::ifstream param_file(param_file_name);
	std::cerr<<param_file_name<<std::endl;	
	if ( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.){
		for(int i = 0 ; i < 6 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfidft1l[i][j];
				param_file >> fgPar_Pfidft1r[i][j];
				param_file >> fgPar_Pfidft2l[i][j];
				param_file >> fgPar_Pfidft2r[i][j];
				param_file >> fgPar_Pfidbt1l[i][j];
				param_file >> fgPar_Pfidbt1r[i][j];
				param_file >> fgPar_Pfidbt2l[i][j];
				param_file >> fgPar_Pfidbt2r[i][j];
				param_file >> fgPar_Pfidbl  [i][j];
				param_file >> fgPar_Pfidbr  [i][j];
			}
		}
    for(int i = 0 ; i < 2 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS2[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS3[i][j];
			}
		}
		for(int i = 0 ; i < 4 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS4[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS5[i][j];
			}
		}
    for(int i = 0 ; i < 2 ; i++){
			for(int j = 0 ; j < 4 ; j++){
				param_file >> fgPar_Pfid_ScpdS2_extra[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 4 ; j++){
				param_file >> fgPar_Pfid_ScpdS3_extra[i][j];
			}
		}
		for(int i = 0 ; i < 4 ; i++){
			for(int j = 0 ; j < 4 ; j++){
				param_file >> fgPar_Pfid_ScpdS4_extra[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 4 ; j++){
				param_file >> fgPar_Pfid_ScpdS5_extra[i][j];
			}
		}
	}
	else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.){
		for(int i = 0 ; i < 6 ; i++){
			for(int j = 0 ; j < 4 ; j++){
				for(int k = 0 ; k < 7 ; k++){
					param_file >> fgPar_Pfid_For[i][j][k];
				}
			}
		}
		for(int i = 0 ; i < 6 ; i++){
			for(int j = 0 ; j < 4 ; j++){
        for(int k = 0 ; k < 7 ; k++){
					param_file >> fgPar_Pfid_Bak[i][j][k];
				}
			}
		}

		for(int i = 0 ; i < 2 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS2[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS3[i][j];
			}
		}
		for(int i = 0 ; i < 4 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS4[i][j];
			}
		}
		for(int i = 0 ; i < 8 ; i++){
			for(int j = 0 ; j < 6 ; j++){
				param_file >> fgPar_Pfid_ScpdS5[i][j];
			}
		}
	}
  else if (E1 == 1161 && torus_current < 1510 && torus_current > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file >> fgPar_1gev_1500_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_1500_Pfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_1500_Pfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_1500_Pfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_1500_Pfid_ScpdS5[i][j];
        }
      }
    }
  else if (E1 == 1161 && torus_current < 760 && torus_current > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file >> fgPar_1gev_750_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_750_Pfid_ScpdS2[i][j];
	  std::cout << fgPar_1gev_750_Pfid_ScpdS2[i][j] << std::endl;
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_750_Pfid_ScpdS3[i][j];
	  std::cout << fgPar_1gev_750_Pfid_ScpdS3[i][j] << std::endl;
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_750_Pfid_ScpdS4[i][j];
	  std::cout << fgPar_1gev_750_Pfid_ScpdS4[i][j] << std::endl;
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file >> fgPar_1gev_750_Pfid_ScpdS5[i][j];
	  std::cout << fgPar_1gev_750_Pfid_ScpdS5[i][j] << std::endl;
        }
      }
    }
	else {
		std::cerr << "read_p_fid_params doesn't have correction parameters for the given input. Check it and fix it!\n";
    exit(-3);
	}
	param_file.close();
	return true;
}
// ===================================================================================================================================

bool Fiducial::read_pim_fid_params()
{
  char param_file_name[256];
  sprintf(param_file_name,"%s/PIMFID_%d_%d.dat",e2adir.c_str(),E1,torus_current);
  std::ifstream param_file(param_file_name);
  if (! param_file.is_open())
    {
      std::cerr << "Failed to open pi- fiducial cut file: \n"
		<< param_file_name << "\n\n Exiting...\n\n";
      exit(-3);
    }
  
  // Read in the pi- parameters

  // Momentum parameters
  for(Int_t sector=0;sector<6;sector++)
    for(Int_t thetapar=0;thetapar<5;thetapar++)
      for(Int_t mompar=0;mompar<6;mompar++)
	param_file >> fgPar_Pimfid[sector][thetapar][mompar];
  
  // S3 parameters
  for(int i = 0 ; i < 4 ; i++)
    for(int j = 0 ; j < 8 ; j++)
      param_file >> fgPar_Pimfid_Theta_S3[i][j];
  
  // S4 Parameters
  for(int i = 0 ; i < 2 ; i++)
    for(int j = 0 ; j < 8 ; j++)
      param_file >> fgPar_Pimfid_Theta_S4[i][j];
  
  // S5 Parameters
  for(int i = 0 ; i < 8 ; i++)
    for(int j = 0 ; j < 8 ; j++)
      param_file >> fgPar_Pimfid_Theta_S5[i][j];
  
  // S3 extra parameters
  for(int i = 0 ; i < 4 ; i++)
    for(int j = 0 ; j < 4 ; j++)
      param_file >> fgPar_Pimfid_Theta_S3_extra[i][j];
  
  // S4 extra parameters
  for(int i = 0 ; i < 2 ; i++)
    for(int j = 0 ; j < 4 ; j++)
      param_file >> fgPar_Pimfid_Theta_S4_extra[i][j];
  
  // S5 extra parameters
  for(int i = 0 ; i < 8 ; i++)
    for(int j = 0 ; j < 4 ; j++)
      param_file >> fgPar_Pimfid_Theta_S5_extra[i][j];

  // We ought to test that the file actually contained all the parameters it needed, but for now...
  return true;
}

bool Fiducial::read_e_fid_params()
{
  char param_file_name[256];
  sprintf(param_file_name,"%s/FCP_%d_%d.dat",e2adir.c_str(),E1,torus_current);
  std::ifstream param_file(param_file_name);
  if(!param_file.is_open()){
    std::cerr<<"Could not open file: "<<param_file_name<<std::endl;
    exit(-1);
  }
  if (E1==4461){
    int param_type=-1;
    int sector=-1;
    double data[6];
    while ( (sector!=6 || param_type!=21))
      {
	param_file >> param_type;
	param_file >> sector >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
	// Test the type of parameter and assign it to the proper data array
	std::cout << param_type << " " << sector << std::endl;
	switch (param_type)
	  {
	  case  0:
	    for(int k=0; k<2; k++) fgPar_Efid_t0_p[sector-1][k] = data[k];
	    break;
	  case  1:
	    for(int k=0; k<6; k++) fgPar_Efid_t1_p[sector-1][k] = data[k];
	    break;
	  case 10:
	    for(int k=0; k<6; k++) fgPar_Efid_b_p[sector-1][0][k] = data[k];
	    break;
	  case 11:
	    for(int k=0; k<6; k++) fgPar_Efid_b_p[sector-1][1][k] = data[k];
	    break;
	  case 20:
	    for(int k=0; k<6; k++) fgPar_Efid_a_p[sector-1][0][k] = data[k];
	    break;
	  case 21:
	    for(int k=0; k<6; k++) fgPar_Efid_a_p[sector-1][1][k] = data[k];
	    break;
	  default:
	    printf("Error in Efid parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	    exit(-1);
	    break;
	  }
      } // Done reading in Fiducial Region Parameters
    // ---
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S3[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S4[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S5[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_Efid_Theta_S3_extra[i][j];
	std::cout <<	fgPar_Efid_Theta_S3_extra[i][j] << std::endl;
      }
    }
    // ---
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_Efid_Theta_S4_extra[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_Efid_Theta_S5_extra[i][j];
      }
    }
  }
  // -----------------------------------------------------
  else if (E1==2261){
    for(int i = 0 ; i < 6 ; i++){
      for(int j = 0 ; j < 6 ; j++){
	for(int k = 0 ; k < 9 ; k++){
	  param_file >> fgPar_Efid[i][j][k];
	}
      }
    }
    // ---
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S3[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S4[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_Efid_Theta_S5[i][j];
      }
    }
  }
  else if (E1 == 1161 && torus_current < 1510 && torus_current > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file >> fgPar_1gev_1500_Efid[sector][thetapar][mompar];
                }
            }
        }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_1500_Efid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_1500_Efid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_1500_Efid_Theta_S5[i][j];
        }
      }
    }
  else if (E1 == 1161 && torus_current < 760 && torus_current > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file >> fgPar_1gev_750_Efid[sector][thetapar][mompar];
                }
            }
        }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_750_Efid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_750_Efid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file >> fgPar_1gev_750_Efid_Theta_S5[i][j];
        }
      }
    }
  // -----------------------------------------------------
  else{
    std::cerr << "File " << param_file_name << " called by Fiducial::read_e_fid_params() does not exist. Check it and fix it!\n";
    exit(-2);
  }
  
  param_file.close();
  
  return true;
}

// ===================================================================================================================================
bool Fiducial::read_vz_cor_params()
{
  // Create file name
  char param_file_name[256];
  sprintf(param_file_name,"%s/vz_%d_%s.dat",e2adir.c_str(),E1,tar.c_str());	

  // open file
  ifstream cal_file(param_file_name);
  
  // test that it's actually open
  if (! cal_file.is_open())
    {
      std::cerr << "File " << param_file_name << " called by Fiducial::read_vz_cor_params() does not exist. Check it and fix it!\n";
      exit(-2);
    }
  
  // Read in the two parameters
  cal_file >> vz_params[0] >> vz_params[1];

  // close the file. 
  cal_file.close();
  
  return true;
}

// ===================================================================================================================================
bool Fiducial::read_e_pcor_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/EMCP_%d_%d.dat",e2adir.c_str(),E1,torus_current);
	std::ifstream param_file(param_file_name);

	int param_type, sector;
	double data[6];
	int cj;
  if (E1 > 2000 && E1 < 5000)
    {
      while (param_file >> param_type)
        {
          param_file >> sector >> cj >> data[0];

          // Assign the data to the correct arrays
          switch (param_type)
            {
            case 0:
              fgPar_Phi[sector-1][cj] = data[0];
              break;
            case 1:
              fgPar_Theta[sector-1][cj] = data[0];
              break;
            default: 
              printf("Error in EMCP parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
              exit(-2);
            }
        }
    }
  else
    {
      for(Int_t eachsecpar=0;eachsecpar<6;eachsecpar++)
        {
          for(Int_t eachphipar=0;eachphipar<6;eachphipar++)
            {
              for(Int_t eachthetapar=0;eachthetapar<4;eachthetapar++)
                {
                  param_file >> fgPar_1gev[eachsecpar][eachphipar][eachthetapar];
                }
            }
        }
    }
    param_file.close();

	// NOTE: corrections for electron momentum obtained with e1c 2.5Gev 2250A data set (Run 16719 and 16720)
	return true;
}
// ===================================================================================================================================
bool Fiducial::read_e_pid_params()
{
	char param_file_name[256];

	if(is_data) sprintf(param_file_name,"%s/el_Epratio_mom_%d.root"    ,e2adir.c_str(),E1);
	else        sprintf(param_file_name,"%s/el_Epratio_mom_%d_sim.root",e2adir.c_str(),E1);

	TFile * old_gfile = gFile;
	TFile * cal_file = new TFile(param_file_name);

	// If we previously set these, we should clean up their memory
	if (el_Ep_ratio_mean)
		delete el_Ep_ratio_mean;
	if (el_Ep_ratio_sig)
		delete el_Ep_ratio_sig;

	// Pull from file
	el_Ep_ratio_mean=(TF1*)cal_file->Get("f_mean")->Clone();
	el_Ep_ratio_sig=(TF1*)cal_file->Get("f_sig")->Clone();

	// Put the root global file pointer back to where it was. I hate ROOT. 
	cal_file->Close();
	gFile = old_gfile;

	// Test that the histograms were pulled successfully
	if (!el_Ep_ratio_mean)
		return false;
	if (!el_Ep_ratio_sig)
		return false;

	return true;
}
// ===================================================================================================================================
bool Fiducial::read_p_pid_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/protdeltat_mom_%d_%d.root",e2adir.c_str(),E1,torus_current);
	TFile * old_gfile = gFile;
	TFile * file_in1 = new TFile(param_file_name);

	// If we previously set these, we should clean up their memory
	if (prot_deltat_sig)	delete prot_deltat_sig;
	if (prot_deltat_mean)	delete prot_deltat_mean;

	// Pull from file
	prot_deltat_sig =(TF1*)file_in1->Get("sig_pol9" )->Clone();
	prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9")->Clone();

	// Put the root global file pointer back to where it was. I hate ROOT. 
	file_in1->Close();
	gFile = old_gfile;

	// Test that the histograms were pulled successfully
	if (!prot_deltat_sig)	return false;
	if (!prot_deltat_mean)	return false;

	return true;
}
// ===================================================================================================================================
bool Fiducial::read_pip_pid_params()
{
        char param_file_name[256];
        sprintf(param_file_name,"%s/pipdeltat_mom_%d_%d.root",e2adir.c_str(),E1,torus_current);
        TFile * old_gfile = gFile;
        TFile * file_in1 = new TFile(param_file_name);

        // If we previously set these, we should clean up their memory
        if (pip_deltat_sig)    delete pip_deltat_sig;
        if (pip_deltat_mean)   delete pip_deltat_mean;

        // Pull from file
        pip_deltat_sig =(TF1*)file_in1->Get("sig_pol9" )->Clone();
        pip_deltat_mean=(TF1*)file_in1->Get("mean_pol9")->Clone();

        // Put the root global file pointer back to where it was. I hate ROOT. 
        file_in1->Close();
        gFile = old_gfile;

        // Test that the histograms were pulled successfully
        if (!pip_deltat_sig)   return false;
        if (!pip_deltat_mean)  return false;

	// 4.4 GeV parameters obtained from runs 179-08,08,10,12,13,14,15,16,17,19,20,21,22
	// 2.2 GeV parameters obtained from runs 181-80,81,82,83,85,86,88,90,98,99 and 182-00,01,02,03,04,05,06

        return true;
}
// ===================================================================================================================================
bool Fiducial::read_pim_pid_params()
{
        char param_file_name[256];
        sprintf(param_file_name,"%s/pimdeltat_mom_%d_%d.root",e2adir.c_str(),E1,torus_current);
        TFile * old_gfile = gFile;
        TFile * file_in1 = new TFile(param_file_name);

        // If we previously set these, we should clean up their memory
        if (pim_deltat_sig)    delete pim_deltat_sig;
        if (pim_deltat_mean)   delete pim_deltat_mean;

        // Pull from file
        pim_deltat_sig =(TF1*)file_in1->Get("sig_pol9" )->Clone();
        pim_deltat_mean=(TF1*)file_in1->Get("mean_pol9")->Clone();

        // Put the root global file pointer back to where it was. I hate ROOT. 
        file_in1->Close();
        gFile = old_gfile;

        // Test that the histograms were pulled successfully
        if (!pim_deltat_sig)   return false;
        if (!pim_deltat_mean)  return false;

        // 4.4 GeV parameters obtained from runs 179-08,08,10,12,13,14,15,16,17,19,20,21,22
        // 2.2 GeV parameters obtained from runs 181-80,81,82,83,85,86,88,90,98,99 and 182-00,01,02,03,04,05,06

        return true;
}
// ===================================================================================================================================
bool Fiducial::in_p_deltaT(double delta_t, double mom, double cut_sigma)
{
	double prot_mom_lim;
	double prot_min = 0.3; //GeV

	if      ( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.) prot_mom_lim=2.;
	else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.) prot_mom_lim=2.;
  else prot_mom_lim=1.;

  if (mom > prot_mom_lim) mom = prot_mom_lim;
	if (mom < prot_min    ) return false;

	double delta_t_up_limit = prot_deltat_mean->Eval(mom) + cut_sigma * prot_deltat_sig->Eval(mom);
	double delta_t_lo_limit = prot_deltat_mean->Eval(mom) - cut_sigma * prot_deltat_sig->Eval(mom);

	if ((delta_t > delta_t_lo_limit) && (delta_t < delta_t_up_limit))
		return true;
	else
		return false;
}
// ===================================================================================================================================
bool Fiducial::in_pip_deltaT(double delta_t, double mom, double cut_sigma)
{
  double pip_mom_lim=2.5;
  
	if      ( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.) pip_mom_lim=2.5;
  else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.) pip_mom_lim=1.4;
  
  if (mom > pip_mom_lim) mom = pip_mom_lim;
  
  double delta_t_up_limit = pip_deltat_mean->Eval(mom) + cut_sigma * pip_deltat_sig->Eval(mom);
  double delta_t_lo_limit = pip_deltat_mean->Eval(mom) - cut_sigma * pip_deltat_sig->Eval(mom);

  if ((delta_t > delta_t_lo_limit) && (delta_t < delta_t_up_limit))
    return true;
  else
    return false;
}
// ===================================================================================================================================
bool Fiducial::in_pim_deltaT(double delta_t, double mom, double cut_sigma)
{
        double pim_mom_lim=2.5;

        if      ( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.) pim_mom_lim=2.5;
        else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.) pim_mom_lim=1.4;

        if (mom > pim_mom_lim) mom = pim_mom_lim;

        double delta_t_up_limit = pim_deltat_mean->Eval(mom) + cut_sigma * pim_deltat_sig->Eval(mom);
        double delta_t_lo_limit = pim_deltat_mean->Eval(mom) - cut_sigma * pim_deltat_sig->Eval(mom);

        if ((delta_t > delta_t_lo_limit) && (delta_t < delta_t_up_limit))
                return true;
        else
                return false;
}
// ===================================================================================================================================
bool Fiducial::in_e_EoverP(double EoverP, double mom, double cut_sigma)
{
	double min_el_mom;
	double max_el_mom;

	if     (E1==4461){
		min_el_mom = 1.10; //GeV
		max_el_mom = 3.70; //GeV
	}
  else if(E1==2261){
		min_el_mom = 0.55; //GeV
		max_el_mom = 2.10; //GeV
	}

  //change maybe?
	else if(E1==1161){
    min_el_mom = 0.30;
    max_el_mom = 1.10;
    }
  else{
		std::cerr << "Fiducial::in_e_EoverP does not have parameters for this beam energy. Check it and fix it!\n";
		exit(-2);
	}

	if (mom < min_el_mom)	return false;
	if (mom > max_el_mom)	mom = max_el_mom;

	double min_EoverP = el_Ep_ratio_mean->Eval(mom) - cut_sigma * el_Ep_ratio_sig->Eval(mom);
	double max_EoverP = el_Ep_ratio_mean->Eval(mom) + cut_sigma * el_Ep_ratio_sig->Eval(mom);

	if ((EoverP > min_EoverP) && (EoverP < max_EoverP))
		return true;
	else
		return false;
}
// ===================================================================================================================================
TVector3 Fiducial::eMomentumCorrection(const TVector3 V3el) const
{
  // Electron Momentum correction, Pass the electron 3 vector, return corrected 3 vector pointer.
  // Check out "http://nuclear.unh.edu/~maurik/Personal/E2Root/html/TE2AnaTool.html"

  if(E1 == 4461){
    return eMomentumCorrection_4GeV(V3el);
  }
  else if(E1 == 2261){ 
    return eMomentumCorrection_2GeV(V3el);
  }
  else if(E1 == 1161){
    return eMomentumCorrection_1GeV(V3el);
  }
  else {
    std::cerr << "eMomentumCorrection doesn't have correction parameters for the given input. Check it and fix it!\n";
  }
  
return V3el;

}

TVector3 Fiducial::eMomentumCorrection_1GeV(TVector3 V3uncor) const{

    // Electron Momentum correction, Pass the electron 3 vector, return corrected 3 vector pointer.
  // Check out "http://nuclear.unh.edu/~maurik/Personal/E2Root/html/TE2AnaTool.html"
  
  double p   = V3uncor.Mag();
  double phi = 180*V3uncor.Phi()/M_PI; 
  if (phi<-30.) phi += 360;
  double theta =  V3uncor.Theta()*180./M_PI;
  int sectInd = (phi+30.)/60.;

  phi -= 60.*sectInd;
  double phipar[6];
  for(int eachphipar=0;eachphipar<6;eachphipar++) {
    phipar[eachphipar]=0;
    for(int eachthetapar=0;eachthetapar<4;eachthetapar++) {
      phipar[eachphipar]+=fgPar_1gev[sectInd][eachphipar][eachthetapar]*pow(theta,eachthetapar);
    }
  }
  Double_t correction=0;
  
  for(int eachphipar=0;eachphipar<6;eachphipar++) {
    correction+=phipar[eachphipar]*pow(phi,eachphipar);
  }
    p*=correction;

  
  TVector3 V3cor(V3uncor);
  V3cor.SetMag(p);
  
  return V3cor;
 
}


TVector3 Fiducial::eMomentumCorrection_2GeV(const TVector3 V3uncor) const {
  
  // Electron Momentum correction, Pass the electron 3 vector, return corrected 3 vector pointer.
  // Check out "http://nuclear.unh.edu/~maurik/Personal/E2Root/html/TE2AnaTool.html"

  // -----------------------------------------------------------------------
  // Correction for Ebeam = 2.2GeV and 2250A data. For more info see:
  // B. Zhang "Electron momentum correction" https://www.jlab.org/Hall-B/secure/e2/bzh/momcorrection.html, 2003.
  // R. Niyazov "Measurement of Correlated Pair Momentum Distributions on 3He(e,e'pp)n with CLAS. PhD thesis, ODU, 2003.
  
  double p   = V3uncor.Mag();
  double cz  = V3uncor.CosTheta();
  double phi = 180*V3uncor.Phi()/M_PI; 
  if (phi<-30.) phi += 360;
  double theta =  V3uncor.Theta()*180./M_PI;
  int sectInd = (phi+30.)/60.;
  
  phi -= 60.*sectInd;
  p = p*(fgPar_Phi[sectInd][0] + fgPar_Phi[sectInd][1]*phi + fgPar_Phi[sectInd][2]*phi*phi);
  
  if     (cz>0.800 && cz<0.885){p=p*(fgPar_Theta[0][0]*sin(fgPar_Theta[0][1]*(cz+fgPar_Theta[0][2])) + fgPar_Theta[0][3]);}
  else if(cz>0.885 && cz<0.935){p=p*(fgPar_Theta[1][0] + fgPar_Theta[1][1]*cz + fgPar_Theta[1][2]*cz*cz);}
  else if(cz>0.935 && cz<0.970){p=p*(fgPar_Theta[2][0] + fgPar_Theta[2][1]*cz + fgPar_Theta[2][2]*cz*cz);}

  TVector3 V3cor(V3uncor);
  V3cor.SetMag(p);
  
  return V3cor;
  
}
  
TVector3 Fiducial::eMomentumCorrection_4GeV(TVector3 V3uncor) const{
  
  // Electron Momentum correction, Pass the electron 3 vector, return corrected 3 vector pointer.
  // Check out "http://nuclear.unh.edu/~maurik/Personal/E2Root/html/TE2AnaTool.html"
  
  // Correction for Ebeam = 4.4GeV and 2250A data (corrections valid only for theta > 16). For more info see:
  // D. Protopopescu "Electron momentum corrections for CLAS at 4.4 GeV" CLAS-NOTE 2001-008 JLAB, 2001.
  // D. Protopopescu "Measurements of a'_{LT} assymetries in (e,e'p) reactions on 4He and 12C," CLAS Analysis Paper, 2003.
  
  double p   = V3uncor.Mag();
  double phi = 180*V3uncor.Phi()/M_PI; 
  if (phi<-30.) phi += 360;  // Corrections are done in degrees, on [-30:330]
  int sectInd = (phi+30.)/60.;

  double theta =  V3uncor.Theta()*180./M_PI; // also in degrees
  if (theta < 16.)
    return V3uncor;

  p *= ((fgPar_Phi[sectInd][0] + fgPar_Phi[sectInd][1]*phi
	 + fgPar_Phi[sectInd][2]*phi*phi)*(fgPar_Theta[sectInd][3]
					   + (fgPar_Theta[sectInd][2]*theta*theta + fgPar_Theta[sectInd][1]*theta
					      + fgPar_Theta[sectInd][0])/exp(theta)));
  
  TVector3 V3cor(V3uncor);
  V3cor.SetMag(p);
  
  return V3cor;
}


bool Fiducial::pFiducialCut(TVector3 momentum){
	//Positive Hadron Fiducial Cut
	//Check out "http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html"
	Bool_t status = kTRUE;

	Float_t theta = momentum.Theta()*180/M_PI;
	Float_t phi   = momentum.Phi()  *180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();

	// ----------------------------------------------------------------------------------------------------------------
	if (E1>4000 && E1<5000 && torus_current>2240. && torus_current<2260.){

		Float_t parfidl [3];    for(Int_t i=0; i<3; i++){parfidl [i]=0;}
		Float_t parfidr [3];    for(Int_t i=0; i<3; i++){parfidr [i]=0;}
		Float_t parfidbl[2];    for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
		Float_t parfidbr[2];    for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
		Float_t cphil =0;       Float_t cphir =0;
		Float_t phi45l=0;       Float_t phi45r=0;
		Float_t phi60l=0;       Float_t phi60r=0;

		Float_t theta_min = 11;
		Float_t theta_max =140;

		bool Forward=kFALSE;            // defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
		Int_t   thetab    =45;          // this variable defines the edge point for Forward<->Backward regions
		Float_t p1        =0.575;       // last bin momentum for region p<0.6 GeV/c
		if(p<0.2) return false;                // momentum less than 0.2 GeV/c, use 0.2 GeV/c
		if(p>4.4) p=4.4;                // momentum greater than 4.4 GeV/c, use 4.4 GeV/c
		// ----------------------------------------------------------
		//Get parametrized values of theta_max for p<0.6 GeV/c region
		if(p<0.6){theta_max=fgPar_Pfidbl[sector][4]+fgPar_Pfidbl[sector][5]*p;}
		//Get parametrized values of theta_max for p>0.6 GeV/c region
		else{theta_max=fgPar_Pfidbl[sector][4]+fgPar_Pfidbl[sector][5]*p1;}

		//Get the momentum dependent parameters for Forward Region (theta <45 deg)   
		Forward=kTRUE;
		if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
			//parameters for hyperbolic function
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_Pfidft1l[sector][j]+fgPar_Pfidft1l[sector][j+1]/p;
				parfidr[i]=fgPar_Pfidft1r[sector][j]+fgPar_Pfidft1r[sector][j+1]/p;
			}
		}
		else{//forward2 defines  regions of momenta and p>0.6 GeV/c
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_Pfidft2l[sector][j]+fgPar_Pfidft2l[sector][j+1]/p;
				parfidr[i]=fgPar_Pfidft2r[sector][j]+fgPar_Pfidft2r[sector][j+1]/p;
			}
		}
		// ----------------------------------------------------------
		phi45l=parfidl [0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
		phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
		if(theta>thetab){//backward region defined by theta >45 deg. 
			if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
			if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

			//Get the momentum dependent parameters for Backward Region

			Forward=kFALSE;
			if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
				//parameters for quadratic function
				for (Int_t i=0; i<3; i++){
					Int_t j=2*i;
					parfidl[i]=fgPar_Pfidbt1l[sector][j]+fgPar_Pfidbt1l[sector][j+1]/p;
					parfidr[i]=fgPar_Pfidbt1r[sector][j]+fgPar_Pfidbt1r[sector][j+1]/p;
				}
				//these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
				for (Int_t i=0; i<2; i++){
					Int_t j=2*i;
					parfidbl[i]=fgPar_Pfidbl[sector][j]+fgPar_Pfidbl[sector][j+1]/p;
					parfidbr[i]=fgPar_Pfidbr[sector][j]+fgPar_Pfidbr[sector][j+1]/p;
				}
			}
			else{//backward2 defines  regions of momenta p>0.6 GeV/c
				//parameters for quadratic function
				for (Int_t i=0; i<3; i++){
					Int_t j=2*i;
					parfidl[i]=fgPar_Pfidbt2l[sector][j]+fgPar_Pfidbt2l[sector][j+1]/p;
					parfidr[i]=fgPar_Pfidbt2r[sector][j]+fgPar_Pfidbt2r[sector][j+1]/p;
				}
				//these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
				for (Int_t i=0; i<2; i++){
					Int_t j=2*i;
					parfidbl[i]=fgPar_Pfidbl[sector][j]+fgPar_Pfidbl[sector][j+1]/p1;
					parfidbr[i]=fgPar_Pfidbr[sector][j]+fgPar_Pfidbr[sector][j+1]/p1;
				}
			}
		}
		// -------------------------------------------------------------------------------------------------------------------
		if(Forward){//Forward region
			if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.   
			cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
			cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
		}
		// -------------------------------------------------------------------------------------------------------------------
		else{//Backward region
			phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
			phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

			if(theta<60){
				cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
				cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
			}
			Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0; 
			//dr and er are theta_flat and phi_edge parameters for phi>0;  
			dl=parfidbl[0];el=parfidbl[1];
			dr=parfidbr[0];er=parfidbr[1];

			if(theta>45&&theta<60){ //BackwardA region
				//try to match parametrized values from Forward region to Backward region parameters
				if(cphil>phi45l)cphil=phi45l;
				if(cphir<phi45r)cphir=phi45r;
			}
			//BackwardB region & phi<0
			else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant 
			else if(theta>dl&&theta<=theta_max){
				cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line 
			else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
			//BackwardB region & phi>0
			if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant 
			else if(theta>dr&&theta<=theta_max){
				cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line 
			else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
		}//Backward Region
		// -------------------------------------------------------------------------------------------------------------------
		if(phi<0) status=(phi>cphil); //check the constrains 
		else if(phi>=0) {status=(phi<cphir);
		}

		if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min
		if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

		//p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for 
		//some range of momentum, where edge does not look good.
		bool s1s4 =(theta<11.7&&(sector==0||sector==3));
		bool s5   =(theta<12.2&& sector==4);
		bool s6   =(theta<11.4&& sector==5);
		if( p>=0.6 && p<1.5 && (s1s4||s5||s6) ) status=kFALSE;

    bool SCpdcut = true;
    if(status && SCpdcut){ // cut bad scintillator paddles
      if(p < 1.0){
        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        //if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar2[i] = badpar2[i]*mom_scpd + fgPar_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
          Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar3[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar3[i] = badpar3[i]*mom_scpd + fgPar_Pfid_ScpdS3[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
          Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
            double temp_mom_scpd = mom_scpd;
          for (Int_t i=0; i<4; i++){
	    if (i==0 || i==1)
              {
                if (mom_scpd > .65)
                  {
                    mom_scpd = .65;
                  }
              }
            else
              mom_scpd = temp_mom_scpd;
	    //std::cout << "temp_mom_scpd = " << temp_mom_scpd << std::endl;
	    //std::cout << "mom_scpd = " << mom_scpd << std::endl;
              
            badpar4[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar4[i] = badpar4[i]*mom_scpd + fgPar_Pfid_ScpdS4[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<2;ipar++){
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
          Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar5[i] = badpar5[i]*mom_scpd + fgPar_Pfid_ScpdS5[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }
      else{
        int tsector = sector + 1;
        double mom_scpd =p;
        // sector 2 has one bad paddles
        if (tsector == 2){
          float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<2; i++){
            badpar2[i] = 0;
            // calculate the parameters using 1/p
            badpar2[i] = fgPar_Pfid_ScpdS2_extra[i][0] + fgPar_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<1;ipar++)
            status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
        }
        if (tsector == 3){
          float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<8; i++){
            badpar3[i] = 0;
            // calculate the parameters using 1/p
            badpar3[i] = fgPar_Pfid_ScpdS3_extra[i][0] + fgPar_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
        // sector 4 has two bad paddle
        else if (tsector == 4){
          float badpar4[4];     // 2 parameters to determine the position of the theta gap
          for (int i=0; i<4; i++){
            badpar4[i] = 0;
            // calculate the parameters using 1/p
            badpar4[i] = fgPar_Pfid_ScpdS4_extra[i][0] + fgPar_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<2;ipar++)
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        // sector 5 has four bad paddles
        else if (tsector == 5){ 
          Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            // calculate the parameters using 1/p
            badpar5[i] = fgPar_Pfid_ScpdS5_extra[i][0] + fgPar_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(Int_t ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }
    return status;
	}
	// ----------------------------------------------------------------------------------------------------------------
	else if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.){

		bool SCpdcut = true;

		Float_t mom_for = p;              // momentum for forward constraints
		if (mom_for<0.3) return false;   // momentum smaller than 300 MeV/c, use 300 MeV/c
		if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
		Float_t mom_bak = p;              // momentum for backward constraints
		if (mom_bak<0.2) return false;   // momentum smaller than 200 MeV/c, use 200 MeV/c
		if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
		Float_t theta0 = 8.5;
		Float_t phi_lower = -24.0;
		Float_t phi_upper =  24.0;
		Float_t par_for[4], par_bak[4];
		for (Int_t i=0; i<4; i++){
			par_for[i] = 0; par_bak[i] = 0;
			for (Int_t d=6; d>=0; d--){
				par_for[i] = par_for[i]*mom_for + fgPar_Pfid_For[sector][i][d];
				par_bak[i] = par_bak[i]*mom_bak + fgPar_Pfid_Bak[sector][i][d];
			}
		}
		if (phi < 0) {
			Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
			status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
		}
		else {
			Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
			status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
		}                     // now the forward constrains are checked
		if ( status ) {       // now check the backward constrains
		  if(theta>par_bak[0]) status = kFALSE;
		  else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && 
					      (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
		}
		if(status && SCpdcut){ // cut bad scintillator paddles
		  Int_t tsector = sector + 1;
		  Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
		  if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
		  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
		  if(tsector==2){      // sector 2 has one bad paddle
		    Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
		    for (Int_t i=0; i<2; i++){
		      badpar2[i] = 0;
		      for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
			  Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
				  badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
				  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
			  Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
          double temp_mom_scpd = mom_scpd;
				for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            {
              if (mom_scpd > .65)
                {
                  mom_scpd = .65;
                }
            }
          else
            mom_scpd = temp_mom_scpd;
          badpar4[i] = 0;
				  for (Int_t d=5; d>=0; d--){
				    badpar4[i] = badpar4[i]*mom_scpd + fgPar_Pfid_ScpdS4[i][d];
				  }                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
				  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
			  Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
			  for (Int_t i=0; i<8; i++){
			    badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
					  badpar5[i] = badpar5[i]*mom_scpd + fgPar_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
			  }
			  for (Int_t ipar=0;ipar<4;ipar++){
			    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
			  }
			}
    }
	}
	else if ( E1 > 1000 && E1 < 2000 && torus_current < 1510 && torus_current > 1490){
    Double_t phipars[5]={0,0,0,0,0};
    bool status = true;
    bool SCpdcut = true;
    if (p < .15)
      p =.15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
      Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
      Int_t tsector = sector + 1;
      if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c

      if(tsector==2){      // sector 2 has one bad paddle
        Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
        for (Int_t i=0; i<2; i++){
          badpar2[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
          }                // calculate the parameters using pol5
        }
        status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
        Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
        for (Int_t i=0; i<8; i++){
          badpar3[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<4;ipar++){
          status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
      }
      else if(tsector==4){ // sector 4 has two bad paddles
        Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
        double temp_mom_scpd = mom_scpd;
        for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            {
              if (mom_scpd > .65)
                {
                  mom_scpd = .65;
                }
            }
          else
            mom_scpd = temp_mom_scpd;
              
          badpar4[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<2;ipar++){
          status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
      }
      else if(tsector==5){ // sector 5 has four bad paddles
        Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
        for (Int_t i=0; i<8; i++){
          badpar5[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<4;ipar++){
          status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
      return status;
    }
    return (status);
  }
  else if ( E1 > 1000 && E1 < 2000 && torus_current < 760 && torus_current > 740){
    Double_t phipars[5]={0,0,0,0,0};
    bool status = true;
    bool SCpdcut = true;
    if (p < .15)
      p = .15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
			if (mom_scpd<0.15)mom_scpd=0.15; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
      if(tsector==2){      // sector 2 has one bad paddle
        Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
        for (Int_t i=0; i<2; i++){
          badpar2[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS2[i][d];
          }                // calculate the parameters using pol5
        }
        status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
        Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
        for (Int_t i=0; i<8; i++){
          badpar3[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS3[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<4;ipar++){
          status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
      }
      else if(tsector==4){ // sector 4 has two bad paddles
        Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
        double temp_mom_scpd = mom_scpd;
        for (Int_t i=0; i<4; i++){
          badpar4[i] = 0;
          if (i==0 || i==1)
            {
              if (mom_scpd > .65)
                {
                  mom_scpd = .65;
                }
            }
          else
            mom_scpd = temp_mom_scpd;
              
          for (Int_t d=5; d>=0; d--){
            badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS4[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<2;ipar++){
          //std::cout << "gap " << ipar+1 << std::endl;
          //std::cout << "momentum " << mom_scpd << std::endl;
          //std::cout << "lower bound " << badpar4[2*ipar] << std::endl;
          //std::cout << "upper bound " << badpar4[2*ipar+1] << std::endl;
          status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
      }
      else if(tsector==5){ // sector 5 has four bad paddles
        Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
        for (Int_t i=0; i<8; i++){
          badpar5[i] = 0;
          for (Int_t d=5; d>=0; d--){
            badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS5[i][d];
          }                // calculate the parameters using pol5
        }
        for (Int_t ipar=0;ipar<4;ipar++){
          status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }
    return status;
  }
// ----------------------------------------------------------------------------------------------------------------
	else{
		std::cerr << "pFiducialCut doesn't have correction parameters for the given input. Check it and fix it!\n";
		exit(-3);
	}
	// ----------------------------------------------------------------------------------------------------------------
	return status;

}
// ===================================================================================================================================
// V_z correction
double Fiducial::vz_corr(TVector3 T3_mom) 
{
	double theta = T3_mom.Theta();
	double phi   = 180./M_PI*T3_mom.Phi();
	double phi_mod = phi+30.;
	if(phi_mod<0) phi_mod+=360;

	//return (-(vz_corr_func->GetParameter(1)))*cos((phi_mod-(vz_corr_func->GetParameter(2)))*M_PI/180.)/tan(theta); 
	return -vz_params[0]*cos((phi_mod - vz_params[1])*M_PI/180.)/tan(theta);

	//Vertex Correction Parameters:
	//E1 = 2.2GeV: obtained from empty run 18283, for 4He
	//E1 = 4.4GeV: obtained from empty run 18522, for 4He (works fine for 3He)

	// Axel: uh hello, this correction is time-dependent, not target or beam energy dependent! This should
	//       eventually get moved to Run_Dependent.cpp
}
// ===================================================================================================================================
TVector3 Fiducial::FindUVW(TVector3 xyz)
{       
	// get the U V W distance to EC edge for the purpose of geometry cut
	// ported from Stepan's function ec_xyz_duvw. the input is lab coordinates of the EC hit.

	Float_t x = xyz.X();
	Float_t y = xyz.Y();
	Float_t z = xyz.Z();
	Float_t xi,yi,zi,u,v,w;
	Float_t ec_the = 0.4363323;
	Float_t ylow   = -182.974 ;
	Float_t yhi    =  189.956 ;
	Float_t tgrho  =  1.95325 ;
	Float_t sinrho = 0.8901256;
	Float_t cosrho = 0.4557150;

	Float_t phi=xyz.Phi()*180./M_PI;
	if(phi<-30) phi+=360;

	Int_t ec_sect = (phi+30)/60.;
	if(ec_sect<0)ec_sect=0;
	if(ec_sect>5)ec_sect=5;

	Float_t ec_phi = ec_sect*M_PI/3.;
	xi = -x*sin(ec_phi) + y*cos(ec_phi);
	yi = x*cos(ec_the)*cos(ec_phi) + y*cos(ec_the)*sin(ec_phi) - z*sin(ec_the);
	zi = x*sin(ec_the)*cos(ec_phi) + y*sin(ec_the)*sin(ec_phi) + z*cos(ec_the);
	zi -= 510.32; 

	u =  ( yi-ylow)/sinrho;
	v =  (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
	w = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;

	TVector3 uvw(u,v,w);
	return uvw;
}       
// ===================================================================================================================================
bool Fiducial::CutUVW_e(TVector3 ecxyz)
{       
	// Cut the edges of EC according to UVW distance threshold defined by par_EcUVW array.
	// If it passes the cut, return true, if not return false

	//parameters for EC edge cuts
  if (E1 >1000 && E1<2000)
    {
      const Float_t par_EcUVW[6][3] = {{20, 368, 403}, {20, 368, 403}, {20, 368, 403}, {20, 368, 403}, {20, 368, 403}, {20, 368, 403}};
      TVector3 ecuvw = FindUVW(ecxyz);
      Float_t phi=ecxyz.Phi()*180/M_PI;
      if(phi<-30) phi+=360;
      Int_t sector = (phi+30)/60;
      if(sector<0)sector=0;
      if(sector>5) sector=5;
      return (ecuvw.X()>par_EcUVW[sector][0] && ecuvw.Y()<par_EcUVW[sector][1] && ecuvw.Z()<par_EcUVW[sector][2]);
    }
  else
    {
      const Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400}, {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};
    
	TVector3 ecuvw = FindUVW(ecxyz);
	Float_t phi=ecxyz.Phi()*180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = (phi+30)/60;
	if(sector<0)sector=0;
	if(sector>5) sector=5;
	return (ecuvw.X()>par_EcUVW[sector][0] && ecuvw.Y()<par_EcUVW[sector][1] && ecuvw.Z()<par_EcUVW[sector][2]);
    }
}
// ===================================================================================================================================
bool Fiducial::CutUVW(TVector3 ecuvw, double dist)
{

        double angle, x_rot, y_rot, d_intercept, abs_y_rot;
        double X_EC = ecuvw.X();
        double Y_EC = ecuvw.Y();
        double n_phi= atan2(Y_EC,X_EC)*180./3.14159;

        if (n_phi < -30.) n_phi += 360.;
        int sec = (int)(n_phi+30)/60;
        if (sec>5) sec = 5;
        if (sec<0) sec = 0;

        angle = 3.14159/180.*60.*(sec);
        d_intercept = dist/cos(atan(1.73));

        x_rot = X_EC*cos(angle) + Y_EC*sin(angle);
        y_rot = X_EC*sin(angle) - Y_EC*cos(angle);

	if(y_rot>=0) abs_y_rot =     y_rot;
	if(y_rot<0)  abs_y_rot = -1.*y_rot;

        return((x_rot<390.-dist)&&(x_rot>1.73*abs_y_rot+55.+d_intercept));
}
// ===================================================================================================================================
double Fiducial::corrected_path_length( double uncorrected_path_length , double E_in , double E_out  ){

        if      (E_in > 0  && E_out == 0)       return uncorrected_path_length + pl_corr_in  ;
        else if (E_in == 0 && E_out > 0 )       return uncorrected_path_length + pl_corr_out ;
        else if (E_in > 0  && E_out > 0 )       return uncorrected_path_length + pl_corr_both;
        else{   std::cerr << "There's something wrong in corrected_path_length" << std::endl;   exit(3);}

}

double Fiducial::EC_in_cut()
{
  if (E1 == 4461)
    return 0.055;
  if (E1 == 2261)
    return 0.060;
  if (E1 == 1161)
    return -999.;
  return -999.;
}

double Fiducial::el_EC_cut()
{
  if (E1==4461)
    return 0.330;
  if (E1==2261)
    return -999.;
  if (E1==1161)
    return -999.;
  return -999.;
}
