#include "Run_dependent.h"
#include "e2a_constants.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TFile.h"
// ===================================================================================================================================
Run_dependent::Run_dependent(int run_number)
{
	// Initialize the key run settings
	run=run_number;

	homedir = std::string(getenv("HOME"));

	// Read in the various parameters
	bool all_ok=true;

	all_ok &= read_run_table();
	all_ok &= read_p_pcor_params(); // proton momentum corrections

	if (all_ok)
		std::cerr << "Successfully read in the various parameters...\n";
	else
	{
		std::cerr << "Failed to read in the parameter files. Exiting...\n";
		exit(-1);
	}  
}
// ===================================================================================================================================
Run_dependent::~Run_dependent()
{
	// Memory clean up
}
// ===================================================================================================================================
bool Run_dependent::read_run_table()
{
  if (run < 10000)
    {
      std::cerr << "This is a simulated run. Making up numbers!!!\n";
      E1=4461;
      torus=2250;
      mini=5996;
      targ="12C";
    }
  else{
    int this_run;
    // Open up the run_table
    std::ifstream run_table((homedir + std::string("/.e2a/run_table.dat")).c_str());
    while (run_table >> this_run)
      {
	run_table >> E1   ;
	run_table >> torus;
	run_table >> mini ;
	run_table >> targ ;
	if (this_run == run)
	  break;
      } 
    
    if (run != this_run)
      {
	std::cerr << "You've selected a run that can't be found on the run table.\n";
	return false;
      }
  }
  return true;
}

bool Run_dependent::read_p_pcor_params() // Proton correction parameters
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/prot_mom_corr_18338_18438.dat",homedir.c_str());
				        
	std::ifstream param_file(param_file_name);

	for(int i = 0 ; i < 6; i++){
		param_file >> up_parm[i];
		param_file >> down_parm[i];
	}
	param_file.close();

	return true;
}

// ===================================================================================================================================
double Run_dependent::ProtonMomCorrection_He3_4Cell(TVector3 V3Pr, double vertex_p ){
	// Low energy proton momentum correction function
	// to be used with He3 target (4th target cell) (RUN # 18338-18438)
	// Input: Proton momentum 3 vector, and Z coord of proton vertex.
	// Returns the corrected MAGNITUDE of the proton momentum,

  double proton_p     = V3Pr.Mag();

  if ((run >= 18338) && (run <= 18438))
    {
      double theta_p      = V3Pr.Theta()*180./M_PI;
      
      double polinom_up   = (((((up_parm[5]*proton_p+up_parm[4])*proton_p+up_parm[3])
			       *proton_p+up_parm[2])*proton_p+up_parm[1])*proton_p+up_parm[0]);
      
      double polinom_down = (((((down_parm[5]*proton_p+down_parm[4])*proton_p+down_parm[3])
			       *proton_p+down_parm[2])*proton_p+down_parm[1])*proton_p+down_parm[0]);
      
      
      if(polinom_up  <0.) polinom_up   = 0;
      if(polinom_down<0.) polinom_down = 0;
      
      double  p_corr_up   = proton_p + proton_p*polinom_up;
      double  p_corr_down = proton_p + proton_p*polinom_down;
      
      if((theta_p>=70.)) return p_corr_up;
      
      if((theta_p < 30.)||(vertex_p>=(1/20.*theta_p-5/2))||(theta_p<=(-200*proton_p+86))) return p_corr_down;
      
      if((theta_p<=70.)&&(theta_p>=30)&&(theta_p>(20*vertex_p+50)))	return p_corr_up;
      else if(proton_p<0.57)						return p_corr_down;
      else 								return p_corr_up;
    }
  
  return proton_p;
}

int Run_dependent::get_E1()
{
  return E1;
}

int Run_dependent::get_torus()
{
  return torus;
}

int Run_dependent::get_mini()
{
  return mini;
}

std::string Run_dependent::get_targ()
{
  return targ;
}
