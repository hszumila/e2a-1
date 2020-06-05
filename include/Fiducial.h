#ifndef __FIDUCIAL_H__
#define __FIDUCIAL_H__

#include <string>
#include "TVector3.h"

class TF1;

class Fiducial
{
	public:
		Fiducial(int E_beam, int torus_current, int mini_current, std::string target, bool data);
		~Fiducial();

		// Functions to test fiducial and pid bounds
		bool e_inFidRegion(TVector3 mom);
		bool pim_inFidRegion(TVector3 mom);
		bool in_e_EoverP(double EoverP, double mom, double cut_sigma);
		bool in_p_deltaT(double delta_t, double mom, double cut_sigma);
		bool in_pip_deltaT(double delta_t, double mom, double cut_sigma);
		bool in_pim_deltaT(double delta_t, double mom, double cut_sigma);
		bool pFiducialCut(TVector3 momentum);
		bool CutUVW_e(TVector3 ecxyz);
		bool CutUVW  (TVector3 ecuvw, double dist);
		double vz_corr(TVector3 T3_mom);
		double corrected_path_length( double uncorrected_path_length , double E_in , double E_out  );
		TVector3 eMomentumCorrection(const TVector3 V3el) const;
		TVector3 FindUVW(TVector3 xyz);
		double EC_in_cut();
		double el_EC_cut();		  
		
	private:
		int E1;
		int torus_current;
		int mini_current;
		std::string e2adir;
		std::string tar;
		bool is_data;

		// Helper functions
		bool read_e_fid_params     ();
		bool read_e_pcor_params    ();
		bool read_e_pid_params     ();
		bool read_p_fid_params     ();
		bool read_pim_fid_params   ();
		bool read_p_pid_params     ();
		bool read_pip_pid_params   ();
		bool read_pim_pid_params   ();
		bool read_vz_cor_params    ();
		bool read_n_pathlength_corr();

		static int get_sector(const TVector3 &p);  //  This returns the 0-indexed sector! Not the one-indexed!

		//Momentum Correction functions are beam energy dependent. Here are the three corrections functions.
		// These specifically are for the original e2a momentum corrections, not Mariana's revised approach...
		TVector3 eMomentumCorrection_1GeV(const TVector3 V3uncor) const;
		TVector3 eMomentumCorrection_2GeV(const TVector3 V3uncor) const;
		TVector3 eMomentumCorrection_4GeV(const TVector3 V3uncor) const;

		// Parameters ///////////////////////////

		// Fiducial Cuts *********
		
		// Electrons
		double fgPar_Efid_t0_p [6][2];
		double fgPar_Efid_t1_p [6][6];
		double fgPar_Efid_b_p  [6][2][6];
		double fgPar_Efid_a_p  [6][2][6];
		double fgPar_Efid      [6][6][9];
		double fgPar_Efid_Theta_S3[4][8]; // Eventually should be replaced by a vector list of gaps
		double fgPar_Efid_Theta_S4[2][8];
		double fgPar_Efid_Theta_S5[8][8];
		double fgPar_Efid_Theta_S3_extra[4][4];
		double fgPar_Efid_Theta_S4_extra[2][4];
		double fgPar_Efid_Theta_S5_extra[8][4];
		double fgPar_1gev_750_Efid[6][5][6];
		double fgPar_1gev_1500_Efid[6][5][6];
		double fgPar_1gev_1500_Efid_Theta_S3[4][8]; // Eventually should be replaced by vector list of gaps
		double fgPar_1gev_1500_Efid_Theta_S4[2][8];
		double fgPar_1gev_1500_Efid_Theta_S5[8][8];
		double fgPar_1gev_750_Efid_Theta_S3[4][8];
		double fgPar_1gev_750_Efid_Theta_S4[2][8];
		double fgPar_1gev_750_Efid_Theta_S5[8][8];
		
		// Positive Hadrons
		double fgPar_Pfidft1l[6][6];
		double fgPar_Pfidft1r[6][6];
		double fgPar_Pfidft2l[6][6];
		double fgPar_Pfidft2r[6][6];
		double fgPar_Pfidbt1l[6][6];
		double fgPar_Pfidbt1r[6][6];
		double fgPar_Pfidbt2l[6][6];
		double fgPar_Pfidbt2r[6][6];
		double fgPar_Pfidbl  [6][6];
		double fgPar_Pfidbr  [6][6];
		double fgPar_Pfid_For[6][4][7];
		double fgPar_Pfid_Bak[6][4][7];
		double fgPar_Pfid_ScpdS2[2][6];
		double fgPar_Pfid_ScpdS3[8][6];
		double fgPar_Pfid_ScpdS4[4][6];
		double fgPar_Pfid_ScpdS5[8][6];
		double fgPar_Pfid_ScpdS2_extra[2][4];
		double fgPar_Pfid_ScpdS3_extra[8][4];
		double fgPar_Pfid_ScpdS4_extra[4][4];
		double fgPar_Pfid_ScpdS5_extra[8][4];
		double fgPar_1gev_750_Pfid[6][5][6];
		double fgPar_1gev_1500_Pfid[6][5][6];
		double fgPar_1gev_1500_Pfid_ScpdS2[2][6];
		double fgPar_1gev_1500_Pfid_ScpdS3[8][6];
		double fgPar_1gev_1500_Pfid_ScpdS4[4][6];
		double fgPar_1gev_1500_Pfid_ScpdS5[8][6];
		double fgPar_1gev_750_Pfid_ScpdS2[2][6];
		double fgPar_1gev_750_Pfid_ScpdS3[8][6];
		double fgPar_1gev_750_Pfid_ScpdS4[4][6];
		double fgPar_1gev_750_Pfid_ScpdS5[8][6];

		// Pi minus 
		double fgPar_Pimfid[6][5][6];
		double fgPar_Pimfid_Theta_S3[4][8];  // Eventually should be replaced by vector list of gaps
		double fgPar_Pimfid_Theta_S4[2][8];
		double fgPar_Pimfid_Theta_S5[8][8];
		double fgPar_Pimfid_Theta_S3_extra[4][4];
		double fgPar_Pimfid_Theta_S4_extra[2][4];
		double fgPar_Pimfid_Theta_S5_extra[8][4];

		// Vertex Corrections ////////////////////////////
		double vz_params[2];

		// Momentum Corrections //////////////////////////
		
		// Electrons **********
		double fgPar_Phi  [6][3];
		double fgPar_Theta[6][4];
		double fgPar_1gev[6][6][4];
		double el_mom_corr_params[6];

		// Particle ID //////////////////////////////////

		// Electrons **********
		TF1 *el_Ep_ratio_mean;
		TF1 *el_Ep_ratio_sig;

		// Protons ************
		TF1 *prot_deltat_sig;
		TF1 *prot_deltat_mean;
		
		// Pi+ ****************
		TF1 *pip_deltat_sig;
                TF1 *pip_deltat_mean;

		// Pi- ****************
                TF1 *pim_deltat_sig;
                TF1 *pim_deltat_mean;

		// Neutron pathlength correction ///////////////
		double pl_corr_in  ;
		double pl_corr_out ;
		double pl_corr_both;
};

#endif
