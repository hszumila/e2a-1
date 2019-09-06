#ifndef __ACCEPTANCE_H__
#define __ACCEPTANCE_H__

#include "TVector3.h"
#include <string>

class TH3D;

class Acceptance
{
 public:
  Acceptance(std::string target, int Ebeam, int torus, std::string part);
  ~Acceptance();
  double get_acc(TVector3 mom);
  double get_acc(double mom, double cosTheta, double phi);
  void get_gen_acc(double mom, double cosTheta, double phi, double &Ngen, double &Nacc);
  void get_bayesian_error_bounds(TVector3 mom, double conf, double &lower, double &upper);

 private:
  TH3D * gen;
  TH3D * acc;
};

#endif
