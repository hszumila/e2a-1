#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>

#include"TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>
#include "TSystem.h"

#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;


int main(int argc, char ** argv)
{
  gSystem->Load("libTree");
	if (argc < 5)
    {
      cerr << "Wrong number of arguments. Instead try\n"
           << "\tcount /path/to/input/file /path/to/output/root /path/to/output/text n or p (looking for which particle)\n\n";
    }

  int no_cut = 0;
  int part_cut = 0;
  int Xb_cut = 0;
  int miss_mom_lower_cut = 0;
  int miss_mom_upper_cut = 0;
  int theta_pq_cut = 0;
  int pq_lower_cut = 0;
  int pq_upper_cut = 0;
  int p_fiducial_cut = 0;
  char type = *argv[argc-1];

  if (type != 'n' && type != 'p')
    {
      cerr << "Wrong type of file. Try instead:\n"
           << "\t (n or p looking for which particle)\n";
      return -2;
    }

  int numfiles = argc-4;
  TFile * infile[numfiles];
  TH1D * vertex[numfiles];
  double mid[numfiles];
  ofstream outtext(argv[numfiles+2]);
  TFile * outfile = new TFile(argv[numfiles+1],"RECREATE");
  for (int file = 0;file<numfiles;file++)
    {
      if(file ==44)
        continue;
      // ---------------------------------------
      // Setting up output tree and branches
      infile[file] = new TFile(argv[file+1]);
      TTree * intree = (TTree*)infile[file]->Get("T");

      int nevents = intree->GetEntries();

      char temp[256];
      sprintf(temp,"file_%d",file+1);
      vertex[file] = new TH1D(temp,temp,100,4,7);

      double e_vz, e_vz_corrected, e_mom[3], e_phi_mod;
      double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
      double EC_in_cut, el_EC_cut;
      double e_t0,beta_assuming_proton,p_t0,delta_t,beta_assuming_pip,pip_t0,pip_delta_t;
      double corr_px, corr_py, corr_pz, n_px, n_py, n_pz, n_p, EC_Path_corr, Beta_corr;

      TVector3 e_ec_xyz, n_ec_xyz;
      TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom, u1;

      int maxPart = 50;
      int nRun, nParticles;
      int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
      int Part_type    [maxPart];
      double Nu, Q2, Xb, Nu_unc, Q2_unc, Xb_unc, t0;
      double vtx_z_unc [maxPart], vtx_z_cor[maxPart], Mass[maxPart];
      double mom_x     [maxPart], mom_y    [maxPart], mom_z  [maxPart]; 
      double e_deltat  [maxPart];
      int    stat_sc   [maxPart], stat_ec  [maxPart], stat_dc[maxPart];
      double sc_time   [maxPart], sc_path  [maxPart];
      double ec_time   [maxPart], ec_path  [maxPart];
      double ec_in     [maxPart], ec_out   [maxPart], ec_tot [maxPart];
      double ec_x      [maxPart], ec_y     [maxPart], ec_z   [maxPart];
      double ec_u      [maxPart], ec_v     [maxPart], ec_w   [maxPart];
      double charge    [maxPart], beta     [maxPart];
      /* 
         =========================
         Part_type values	(http://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html)
         -11  = electron
         2212 = proton
         2112 = neutron
         +211 = pi+
         -211 = pi-
         +111 = pi0
         =========================
      */ 

      intree->SetBranchAddress("nRun"      , &nRun      );
      intree->SetBranchAddress("nParticles", &nParticles);
      intree->SetBranchAddress("nProtons"  , &nProtons  );
      intree->SetBranchAddress("nNeutrons" , &nNeutrons );
      //intree->SetBranchAddress("nPiplus"   , &nPiplus   );
      //intree->SetBranchAddress("nPiminus"  , &nPiminus  );
      //intree->SetBranchAddress("t0"        , &t0        );
      intree->SetBranchAddress("Nu"        , &Nu        );
      //intree->SetBranchAddress("Q2"        , &Q2        );
      //intree->SetBranchAddress("Xb"        , &Xb        );
      //intree->SetBranchAddress("charge"    ,  charge    );
      //intree->SetBranchAddress("beta"      ,  beta      );
      intree->SetBranchAddress("Part_type" ,  Part_type );
      //intree->SetBranchAddress("vtx_z_unc" ,  vtx_z_unc );
      intree->SetBranchAddress("vtx_z_cor" ,  vtx_z_cor );
      intree->SetBranchAddress("mom_x"     ,  mom_x     );
      intree->SetBranchAddress("mom_y"     ,  mom_y     );
      intree->SetBranchAddress("mom_z"     ,  mom_z     );
      /*intree->SetBranchAddress("e_deltat"  ,  e_deltat  );
      intree->SetBranchAddress("stat_sc"   ,  stat_sc   );
      intree->SetBranchAddress("stat_dc"   ,  stat_dc   );
      intree->SetBranchAddress("stat_ec"   ,  stat_ec   );
      intree->SetBranchAddress("sc_time"   ,  sc_time   );
      intree->SetBranchAddress("sc_path"   ,  sc_path   );
      intree->SetBranchAddress("ec_time"   ,  ec_time   );
      intree->SetBranchAddress("ec_path"   ,  ec_path   );
      intree->SetBranchAddress("ec_in"     ,  ec_in     );
      intree->SetBranchAddress("ec_out"    ,  ec_out    );
      intree->SetBranchAddress("ec_tot"    ,  ec_tot    );
      intree->SetBranchAddress("ec_x"      ,  ec_x      );
      intree->SetBranchAddress("ec_y"      ,  ec_y      );
      intree->SetBranchAddress("ec_z"      ,  ec_z      );
      intree->SetBranchAddress("ec_u"      ,  ec_u      );
      intree->SetBranchAddress("ec_v"      ,  ec_v      );
      intree->SetBranchAddress("ec_w"      ,  ec_w      );
      intree->SetBranchAddress("Mass"      ,  Mass      );
      */// --------------------------------------------------------------------------------------------------

      intree->GetEvent(0);
      int tab_run, tab_E1, tab_torus, tab_mini;
      string tab_targ;
      char param_file_name[256];
      string homedir = string(getenv("HOME"));
      sprintf(param_file_name,"%s/.e2a/run_table.dat",homedir.c_str());
      ifstream run_table;
      run_table.open(param_file_name);

      do{
        run_table >> tab_run  ;
        run_table >> tab_E1   ;
        run_table >> tab_torus;
        run_table >> tab_mini ;
        run_table >> tab_targ ;
      } while(tab_run != nRun);

      cout << "Run    = " << tab_run   << endl;
      cout << "Ebeam  = " << tab_E1    << endl;
      cout << "Torus  = " << tab_torus << endl;
      cout << "Mini   = " << tab_mini  << endl;
      cout << "Target = " << tab_targ  << endl;


      int elec_index = 0;
      double beam_energy = tab_E1/1000.;
      TVector3 q,elec_mom,miss_mom,part_mom;

      cout << "Before fiducial class" << endl;
      Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, true);
      // Create an instance of the Fiducial Class

      cout << "After fiducial class" << endl;

      for (int event = 0; event < nevents; event++)
        {
          bool no_cut_bool = false;
          bool part_cut_bool = false;
          bool Xb_cut_bool = false;
          bool miss_mom_lower_cut_bool = false;
          bool miss_mom_upper_cut_bool = false;
          bool theta_pq_cut_bool = false;
          bool pq_lower_cut_bool = false;
          bool pq_upper_cut_bool = false;
          bool p_fiducial_cut_bool = false;
          //if (event%100000==0)
          //cout << "Running through event " << event << endl;

          intree->GetEntry(event);

          TVector3 elec_mom(mom_x[elec_index],mom_y[elec_index],mom_z[elec_index]);
          TVector3 q = (TVector3(0,0,beam_energy) - elec_mom);
          double omega = beam_energy - elec_mom.Mag();
          double q_squared = q.Mag2()-omega*omega;
          double xb = q_squared/(2*0.938272*omega);


          for (int part = 0; part<nParticles;part++)
            {
              if (type == 'n')
                {
                  if (Part_type[part] != 2112)
                    {
                      continue;
                    }
                }

              if (type == 'p')
                {
                  if (Part_type[part] != 2212)
                    {
                      continue;
                    }
                }
              vertex[file]->Fill(vtx_z_cor[part]);
              TVector3 part_mom(mom_x[part],mom_y[part],mom_z[part]);
              TVector3 miss_mom = part_mom-q;
              double miss_mom_mag = miss_mom.Mag();

              part_cut_bool = true;

              if (xb < 1.2)
                continue;

              Xb_cut_bool = true;

              double theta_pq = q.Angle(part_mom)*(180./3.141592);

              if (theta_pq > 25.)
                continue;

              theta_pq_cut_bool = true;


              if(part_mom.Mag()/q.Mag() < .62)
                continue;

              pq_lower_cut_bool = true;

              if(part_mom.Mag()/q.Mag() > .96)
                continue;

              pq_upper_cut_bool = true;


              if (miss_mom_mag < .3)
                continue;

              miss_mom_lower_cut_bool = true;

              if (miss_mom_mag > 1.0)
                continue;

              miss_mom_upper_cut_bool = true;

              // Passing positive hadron fiducial cuts
              if(!fid_params.pFiducialCut(part_mom))
                continue;

              p_fiducial_cut_bool = true;
              /*cout << "The particle is a " << Part_type[part] << endl;
              cout << "Xb is " << xb << endl;
              cout << "The missing momentum is " << miss_mom_mag << endl;
              cout << "The angle between the particle and q is " << theta_pq << endl;8*/
            }
          if (part_cut_bool == true)
            part_cut++;

          if (Xb_cut_bool == true)
            Xb_cut++;

          if (miss_mom_lower_cut_bool == true)
            miss_mom_lower_cut++;

          if (miss_mom_upper_cut_bool == true)
            miss_mom_upper_cut++;

          if (theta_pq_cut_bool == true)
            theta_pq_cut++;

          if (pq_lower_cut_bool == true)
            pq_lower_cut++;

          if (pq_upper_cut_bool == true)
            pq_upper_cut++;

          if (p_fiducial_cut_bool == true)
            p_fiducial_cut++;
        }
      cout << " " << endl;
      cout << "The total number of " << type << " is " << part_cut << endl;
      cout << "The total number passing Xb cuts is " << Xb_cut << endl;
      cout << "after thetapq " << theta_pq_cut << endl;
      cout << "after p/q lower " << pq_lower_cut << endl;
      cout << "after p/q upper " << pq_upper_cut << endl;
      cout << "The total number passing lower missing momentum cuts is " << miss_mom_lower_cut << endl;
      cout << "The total number passing upper missing momentum cuts is " << miss_mom_upper_cut << endl;
      cout << "The total number passing proton fiducial cuts " << p_fiducial_cut << endl;
      TF1 * fit1 = new TF1("tot","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*exp(-((x-[6])*(x-[6]))/(2*[7]*[7]))", -100, 100);
      fit1->SetParameters(200,50,2,.5,.5,vertex[file]->GetMaximum(),5,.04);

      vertex[file]->Fit("tot","","",4,7);
      mid[file] = vertex[file]->GetFunction("tot")->GetParameter(6);
      outtext << tab_run << " " << mid[file] << endl;
      outfile->cd();
      vertex[file]->Write();
      delete vertex[file];
      infile[file]->Close();
    }

  outtext.close();
  outfile->Close();
}
