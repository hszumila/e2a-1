#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "TVector3.h"


#include "Acceptance.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 5)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
           << "\tmap_query [target] [E_beam (MeV)] [torus (A)] [particle]\n\n"
           << "Query points are passed in by stdin in triplets of\n"
           << "\t[mom (GeV/c)] [theta (deg.)] [phi (deg.)]\n\n";
      return -1;
    }
  
  Acceptance myMap(string(argv[1]),atoi(argv[2]),atoi(argv[3]),string(argv[4]));

  double mom, theta_deg, phi_deg;
  while (cin >> mom)
    {
      cin >> theta_deg >> phi_deg;

      double theta = theta_deg * M_PI/180.;
      double phi = phi_deg * M_PI/180.;
      TVector3 temp(mom*sin(theta)*cos(phi), mom*sin(theta)*sin(phi), mom*cos(theta));
      cout << mom << " " << theta_deg << " " << phi_deg << " " << myMap.get_acc(temp) << "\n";
    }

  return 0;
}
