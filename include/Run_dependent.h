#ifndef __RUN_DEPENDENT_H__
#define __RUN_DEPENDENT_H__

#include <string>
#include "TVector3.h"

class TF1;

class Run_dependent
{
	public:
		Run_dependent(int run_number);
		~Run_dependent();
		int get_E1();
		int get_torus();
		int get_mini();
		std::string get_targ();

		double ProtonMomCorrection_He3_4Cell(TVector3 V3Pr, double vertex_p );

	private:
		// Key information about the run
		int run;
		int E1;
		int torus;
		int mini;
		std::string targ;	
		std::string homedir;

		// Helper functions
		bool read_run_table();
		bool read_p_pcor_params();

		// Proton Momentum Correction Data
		double up_parm[6];
		double down_parm[6];


};

#endif
