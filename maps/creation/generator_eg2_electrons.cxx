#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv){

  if (argc != 2)
    {
      cerr << "generator is a program that produces a file 'mctk_uniform.txt'\n"
	   << "  for use with gsim. You must supply a phi angle (1-360).\n\n"
	   << "  generator [phi]\n\n";
      return -1;
    }
  else
    {
      cerr << "You have called the generator with phi = " << argv[1] << "\n";
    }

  ofstream outfile;
  outfile.open ("./mctk_uniform.txt");

  Int_t run_number = atoi(argv[1]);

  //Number of Events to Generate
  Int_t num_p_bins = 100;
  Int_t num_theta_bins = 200;
  Int_t num_generated = 10;
  Int_t nEntries = num_p_bins*num_theta_bins*num_generated;
  cout << "Total Number of Entries = " << nEntries << endl;
  
  // const Int_t minCosTBin = (Int_t)((1.+TMath::Cos(71.*TMath::DegToRad()))/2. * ((double)num_cost_bins));
  
  //Set number of particles to generate (order:e,p,pi+,pi-)

  // Get seed from /dev/urandom                                                                   
  unsigned int seed;
  FILE *f=fopen("/dev/urandom","r");
  fread(&seed,sizeof(unsigned int),1,f);
  fclose(f);

  //Get Random Number Generator                                                                    
  TRandom3 *myRandom = new TRandom3(seed);

  cout << "The Random Number seed is "<<myRandom->GetSeed()<<endl;

  //Set Variables
  Float_t cx[50],cy[50],cz[50],mom_tot[50];
  Int_t pid[50];Float_t mass[50];Int_t charge[50];

  //These Don't Change
  Int_t pid_e = 11; Float_t mass_e = 0.0005; Int_t charge_e = -1; //electron
  Int_t pid_p = 2212; Float_t mass_p = 0.9383; Int_t charge_p = 1; //proton
  Int_t pid_pp = 211; Float_t mass_pp = 0.1396; Int_t charge_pp = 1; //piplus
  Int_t pid_pm = -211; Float_t mass_pm = 0.1396; Int_t charge_pm = -1; //piminums
  Float_t x = 0.000, y = 0.000, z = -0.000;
  Float_t t_off = 0.000;
  Int_t flag = 0;

  //Set temp variables
  Double_t theta , phi; //i.e. Cos(theta),Phi
  Float_t px, py, pz;

  //Loop over Data
  for(Int_t p_bin=0;p_bin<num_p_bins;p_bin++)
    {
      for(Int_t theta_bin=0;theta_bin<num_theta_bins;theta_bin++)
	{
	  for(Int_t gen=0;gen<num_generated;gen++)
	    {
	      //Get Info for the electron
	      mom_tot[0] = myRandom->Uniform(.05*p_bin, .05*p_bin + .05);
	    
	      theta = myRandom->Uniform((5+.1*theta_bin)*TMath::Pi()/180, (5.1+.1*theta_bin)*TMath::Pi()/180);
	      phi =  myRandom->Uniform(2*TMath::Pi()*(run_number-1)/360., (2*TMath::Pi()*run_number)/360.);

	      px = mom_tot[0] * TMath::Sin(theta) * TMath::Cos( phi );
	      py = mom_tot[0] * TMath::Sin(theta) * TMath::Sin( phi );
	      pz = mom_tot[0] * TMath::Cos(theta); 

	      cx[0] = px/mom_tot[0]; cy[0] = py/mom_tot[0]; cz[0] = pz/mom_tot[0];
	      pid[0] = pid_e; mass[0] = mass_e; charge[0] = charge_e;
	      //------------------------------------------------------------------------
	      for(Int_t k=0;k<1;k++){      
		if(k==0) {outfile << "1" << endl;}
		outfile << pid[k] <<" "<< cx[k] <<" "<< cy[k] <<" "<< cz[k] <<" "<< mom_tot[k] <<endl;
		outfile << mass[k] <<" "<< charge[k] << endl;
		outfile << x <<" "<< y <<" "<< z <<" "<< t_off <<" "<< flag <<endl;
	      }
	    }
	}
    }
  //------------------------------------------------------------------------
  outfile.close();
  
  return 0;
}
  
