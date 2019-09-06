#include "Acceptance.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "TFile.h"
#include "TH3.h"
#include "TEfficiency.h"

Acceptance::Acceptance(std::string target, int Ebeam, int torus, std::string part)
{
  // Sanitize target name
  if (target == "C12") target = "solid";
  if (target == "Fe56") target = "solid";

  // We need to construct the file name
  char * home=getenv("HOME");
  char filename[100];
  sprintf(filename,"%s/.e2a/maps/e2a_%s_%d_%d_%s.root",home, target.c_str(),Ebeam,torus,part.c_str());
  
  // Test if the file exists
  TFile * f = NULL;
  f = new TFile(filename);
  if (f->IsZombie())
    {
      std::cerr << "Attempts to open the file at \n\t" << filename << "\n failed. Double check!\n"
                << "\t It is possible that this map hasn't been created yet. Exiting...\n";
      exit(-1); 
    }

  // Get histograms from file
  gen = (TH3D*)f->Get("Generated Particles");
  acc = (TH3D*)f->Get("Accepted Particles");

  if ((!gen) || (!acc))
    {
      std::cerr << "Couldn't get maps from file " << filename << ".\n"
                << "\tExiting...\n\n";
      exit (-2);
    }

  // Decouple them from the file
  gen->SetDirectory(0);
  acc->SetDirectory(0);

  // Close the open root file.
  f->Close();
}

Acceptance::~Acceptance()
{
  delete gen;
  delete acc;
}

double Acceptance::get_acc(TVector3 mom)
{
  return get_acc(mom.Mag(), mom.CosTheta(), mom.Phi());
}

double Acceptance::get_acc(double mom, double cosTheta, double phi)
{ 
  double Ngen, Nacc;
  get_gen_acc(mom, cosTheta, phi, Ngen, Nacc);

  /*if (Ngen <= 0.)
    {
      std::cerr << "There is an error in this map! The number of generated events in this bin is zero.\n"
		<< "\tPrinting info:\n"
		<< "\t\tMom: " << mom << " theta: " << acos(cosTheta)*180./M_PI << " phi: " << phi*180./M_PI << "\n"
		<< "\t Returning zero acceptance. Please double check.\n";
    return 0.;
    }*/
  
  return Nacc/Ngen;
}

void Acceptance::get_gen_acc(double mom, double cosTheta, double phi, double &Ngen, double &Nacc)
{
  double phi_deg=phi*180./M_PI;

  // Do a range test
  if ((cosTheta < -1.) || (cosTheta > 1.))
    {
      std::cerr << "Provided theta (" << acos(cosTheta)*180./M_PI << " deg.) is out of range. Returning 0 acceptance. Please check.\n";
      Nacc=0.;
      Ngen=10.;
      return;
    }
  if (mom < 0.)
    {
      std::cerr << "Provided momentum (" << mom << " GeV/c) is negative. Returning 0 acceptance. Please check.\n";
      Nacc=0.;
      Ngen=10.;
      return;
    }
  
  // Sanitize momentum and phi
  if (cosTheta == 1.) cosTheta = 0.99999; // Necessary because ROOT hist bins don't include upper bound.
  if (mom >= 5.) mom=4.9999;

  while (phi_deg < -30.)
    phi_deg += 360.;
  while (phi_deg >= 330.)
    phi_deg -= 360.;
  
  int pbin = gen->GetXaxis()->FindBin(mom);
  int cbin = gen->GetYaxis()->FindBin(cosTheta);
  int fbin = gen->GetZaxis()->FindBin(phi_deg);

  Ngen = gen->GetBinContent(pbin,cbin,fbin);
  Nacc = acc->GetBinContent(pbin,cbin,fbin);
}

void Acceptance::get_bayesian_error_bounds(TVector3 mom, double conf, double &lower, double &upper)
{
  double p = mom.Mag();
  double cosTheta = mom.CosTheta();
  double phi = mom.Phi();
  double Ngen, Nacc;
  get_gen_acc(p, cosTheta, phi, Ngen, Nacc);

  if (Ngen <= 0.)
    {
      std::cerr << "There is an error in this map! The number of generated events in this bin is zero.\n"
                << "\tPrinting info:\n"
                << "\t\tMom: " << p << " theta: " << acos(cosTheta)*180./M_PI << " phi: " << phi*180./M_PI << "\n"
                << "\t Returning bounds at 0 and 1. Please double check.\n";
      lower=0.;
      upper=1.;
      return;
    }

 lower=TEfficiency::Bayesian(Ngen, Nacc, conf, 1., 1.,false, true);
 upper=TEfficiency::Bayesian(Ngen, Nacc, conf, 1., 1.,true, true);
}
