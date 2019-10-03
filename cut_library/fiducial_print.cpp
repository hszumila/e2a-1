#include <iostream>
#include <cmath>
#include <cstdlib>
#include "TVector3.h"

#include "Fiducial.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 5)
    {
      cerr << "Wrong number of arguments. Instead use: \n"
	   << "\tfiducial_test [E_beam (MeV)] [torus (A)] [target (ie. 12C)] [sector (0-5)]\n\n";
      return -1;
    }

  int E_beam = atoi(argv[1]);
  if (!((E_beam == 4461)||(E_beam == 2261)||(E_beam == 1161)))
    {
      cerr << "Invalid beam energy provided. Try 4461, 2261 or 1161!\n\n";
      return -2;
    }

  Fiducial my_Fid(E_beam,atoi(argv[2]),5996,string(argv[3]),true);

  const int sec = atoi(argv[4]);
  const double midPhiDeg = 60. * sec;
  const double midPhi = midPhiDeg * M_PI/180.;

  // Loop in momentum
  for (double p=1.25 ; p<5 ; p+=0.25)
    {
      //cerr << "Working on p=" << p << "\n";
      // Find theta min
      double thetaMin=0.;
      while (!my_Fid.e_inFidRegion(TVector3(p*sin(thetaMin)*cos(midPhi),p*sin(thetaMin)*sin(midPhi),p*cos(thetaMin))))
	{
	  thetaMin += 0.0001;
	  //cerr << "Attempting thetaMin = " << thetaMin*180./M_PI <<" with result " << my_Fid.e_inFidRegion(TVector3(p*sin(thetaMin)*cos(midPhi),p*sin(thetaMin)*sin(midPhi),p*cos(thetaMin))) << "\n";
	  if (thetaMin > M_PI/3.)
	    break;
	}

      // Step over theta from theta min to theta max
      for (double theta=thetaMin ; theta < M_PI/3. ; theta+=0.0001)
	{
	  cout << p << " " << theta*180./M_PI << " ";

	  double maxPhiGuess=0.;
	  double minPhiGuess=-M_PI/6.;
	  while (maxPhiGuess - minPhiGuess > 0.0001)
	    {
	      double avgPhi = 0.5 * (maxPhiGuess + minPhiGuess);
	      TVector3 trial(p*sin(theta)*cos(avgPhi + midPhi), p*sin(theta)*sin(avgPhi+midPhi),p*cos(theta));
	      if (my_Fid.e_inFidRegion(trial))
		maxPhiGuess = avgPhi;
	      else
		minPhiGuess = avgPhi;
	    }
	  cout << 0.5*(minPhiGuess + maxPhiGuess)*180./M_PI << " ";

	  maxPhiGuess=M_PI/6.;
	  minPhiGuess=0.;
	  while (maxPhiGuess - minPhiGuess > 0.0001)
	    {
	      double avgPhi = 0.5 * (maxPhiGuess + minPhiGuess);
	      TVector3 trial(p*sin(theta)*cos(avgPhi + midPhi), p*sin(theta)*sin(avgPhi+midPhi),p*cos(theta));
	      if (my_Fid.e_inFidRegion(trial))
		minPhiGuess = avgPhi;
	      else
		maxPhiGuess = avgPhi;
	    }
	  cout << 0.5*(minPhiGuess + maxPhiGuess)*180./M_PI << "\n";
	}
      cout << "\n\n";
    }

  return 0;
}
