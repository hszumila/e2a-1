#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TRint.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv){

  ofstream outfile;
  outfile.open ("./mctk_uniform.txt");

  //Number of Events to Generate
  Int_t nEntries = 100000;
  cout << "Total Number of Entries = " << nEntries << endl;
  
  //Set number of particles to generate (order:e,p,pi+,pi-)
  Int_t top_num = 4;

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
  Float_t x = 0.000, y = 0.000, z = 0.000 ;
  Float_t t_off = 0.000;
  Int_t flag = 0;

  //Set temp variables
  Double_t cost , phi; //i.e. Cos(theta),Phi
  Float_t px, py, pz;

  //Loop over Data
  for(Int_t i=0;i<nEntries;i++){
    
    if(i%100000==0) cout << " events processed = " << i << endl;
    
    //Get Info for the electron
    mom_tot[0] = myRandom->Uniform(0,5);
    
    //To generate uniformly, we should do the following: "cos(theta) = 1 - 2*Uniform[0,1]"
    //This is the same as "cos(theta) = Uniform[-1,1]
    //For electrons restrict to forward hemisphere (0-71 degrees), since no electron detectors in back
    cost = myRandom->Uniform(0.325,1);
    phi =  myRandom->Uniform(0, 2*TMath::Pi());

    px = mom_tot[0] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
    py = mom_tot[0] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
    pz = mom_tot[0] * cost; 

    cx[0] = px/mom_tot[0]; cy[0] = py/mom_tot[0]; cz[0] = pz/mom_tot[0];
    pid[0] = pid_e; mass[0] = mass_e; charge[0] = charge_e;
    
    //Get Info for the proton
    mom_tot[1] = myRandom->Uniform(0,3.5); //Hard cut at reconstructed momentum = 2.8
    cost = myRandom->Uniform(-1,1);
    phi =  myRandom->Uniform(0, 2*TMath::Pi());

    px = mom_tot[1] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
    py = mom_tot[1] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
    pz = mom_tot[1] * cost; 

    cx[1] = px/mom_tot[1]; cy[1] = py/mom_tot[1]; cz[1] = pz/mom_tot[1];
    pid[1] = pid_p; mass[1] = mass_p; charge[1] = charge_p;

    //Get Info for the pi+
    mom_tot[2] = myRandom->Uniform(0,4);
    cost = myRandom->Uniform(-1,1);
    phi =  myRandom->Uniform(0, 2*TMath::Pi());

    px = mom_tot[2] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
    py = mom_tot[2] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
    pz = mom_tot[2] * cost;

    cx[2] = px/mom_tot[2]; cy[2] = py/mom_tot[2]; cz[2] = pz/mom_tot[2];
    pid[2] = pid_pp; mass[2] = mass_pp; charge[2] = charge_pp;

    //Get Info for the pi-
    mom_tot[3] = myRandom->Uniform(0,4);
    cost = myRandom->Uniform(-1,1);
    phi =  myRandom->Uniform(0, 2*TMath::Pi());

    px = mom_tot[3] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
    py = mom_tot[3] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
    pz = mom_tot[3] * cost;

    cx[3] = px/mom_tot[3]; cy[3] = py/mom_tot[3]; cz[3] = pz/mom_tot[3];
    pid[3] = pid_pm; mass[3] = mass_pm; charge[3] = charge_pm;
  
    //------------------------------------------------------------------------
    for(Int_t k=0;k<top_num;k++){
      
      if(k==0) {outfile << top_num << endl;}
      outfile << pid[k] <<" "<< cx[k] <<" "<< cy[k] <<" "<< cz[k] <<" "<< mom_tot[k] <<endl;
      outfile << mass[k] <<" "<< charge[k] << endl;
      outfile << x <<" "<< y <<" "<< z <<" "<< t_off <<" "<< flag <<endl;
    }
    //------------------------------------------------------------------------
    
  }
  
  outfile.close();
  return 0;
}
