#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>

#include"TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>

#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

// Difference between positive corrected z vertex and electron corrected z vertex
const double pos_z_cut_min = -2;
const double pos_z_cut_max =  2;

// Parameters for 2GeV data cuts
const double sc_cc_dt_cut_sect[6]={-2,-5,-8,-8,-2,2};

int main(int argc, char ** argv){
  if (argc < 5){
    cerr << "Wrong number of arguments. Instead try\n"
         << "\tskim_tree /path/to/output/file 00000000 00000 /path/to/input1 [optional: /path/to/input2 ...]\n\n"
         << "Argument 2 must be an 8 digit number:\n"
         << "\t1st digit: minimum number of protons required in the event (0-9)\n"
         << "\t2nd digit: maximum number of protons required in the event (0-9)\n"
         << "\t3rd digit: minimum number of neutrons required in the event (0-9)\n"
         << "\t4th digit: maximum number of neutrons required in the event (0-9)\n"
         << "\t5th digit: minimum number of pi+ required in the event (0-9)\n"
         << "\t6th digit: maximum number of pi+ required in the event (0-9)\n"
         << "\t7th digit: minimum number of pi- required in the event (0-9)\n"
         << "\t8th digit: maximum number of pi- required in the event (0-9)\n\n"
	 << "\tIf any maximum is smaller than the minimum, that particle type will be ignored entirely.\n\n"
         << "Argument 3 must be a 5 digit sequence of 1,0:\n"
         << "\t1st digit: fiducial cut status for electrons (0 = off, 1 = on)\n"
         << "\t2nd digit: fiducial cut status for protons   (0 = off, 1 = on)\n"
         << "\t3rd digit: fiducial cut status for neutrons  (0 = off, 1 = on)\n"
         << "\t4th digit: fiducial cut status for pi+       (0 = off, 1 = on)\n"
         << "\t5th digit: fiducial cut status for pi-       (0 = off, 1 = on)\n\n";
    return -1;
  }

  // --------------------------------------------------------------------------------------------------
  // Specifying the conditions for the output file
  if(strlen(argv[2]) != 8){
    cout << endl << "First argument must be an integer of 8 digits." << endl;
    cout << "You have specified an integer of " << strlen(argv[2]) << " digits" << endl;
    return -2;
  }

  if(strlen(argv[3]) != 5){
    cout << endl << "Second argument must be an integer of 5 digits." << endl;
    cout << "You have specified an integer of " << strlen(argv[3]) << " digits" << endl;
    return -3;
  }

  int in_num_part = atoi(argv[2]); 
  const int min_p   = (in_num_part/10000000)%10;
  const int max_p   = (in_num_part/1000000 )%10;
  const int min_n   = (in_num_part/100000  )%10;
  const int max_n   = (in_num_part/10000   )%10;
  const int min_pip = (in_num_part/1000    )%10;
  const int max_pip = (in_num_part/100     )%10;
  const int min_pim = (in_num_part/10      )%10;
  const int max_pim = (in_num_part/1       )%10;

  int fid_settings = atoi(argv[3]);
  const bool do_fid_e   = fid_settings/10000 & 1;
  const bool do_fid_p   = fid_settings/1000 & 1;
  const bool do_fid_n   = fid_settings/100 & 1;
  const bool do_fid_pip = fid_settings/10 & 1;
  const bool do_fid_pim = fid_settings & 1;

  cout << endl << "The output root file will be skimmed for reactions containing the following number of particles" << endl;
  cout << min_p   << " <= number of protons  <= " << max_p   << endl;
  cout << min_n   << " <= number of neutrons <= " << max_n   << endl;
  cout << min_pip << " <= number of pion+    <= " << max_pip << endl;
  cout << min_pim << " <= number of pion-    <= " << max_pim << endl;	

  // --------------------------------------------------------------------------------------------------
  // Open up the input root files
  TChain t("data");
  for (int i=4 ; i<argc ; i++)
    {
      TFile *temp = new TFile(argv[i]);
      bool file_ok=true;
      if (!temp)
        file_ok = false;
      else if (temp->IsZombie())
        file_ok =false;

      if (file_ok)
        {
          cerr << "Successfully opened file " << argv[i] << "\n";
          temp->Close();
          t.Add(argv[i]);
        }
      else
        {
          cerr << "Could not open file " << argv[i] << "\n"
               << " ... exiting.\n";
          return -4;
        }
    }

  // --------------------------------------------------------------------------------------------------
  // Open up the tree, and get the important data
  const int nEvents = t.GetEntries();
  const int maxPart = 50;
  int gPart, CCPart, DCPart, ECPart, SCPart, NRun;
  int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart], SC_Pad[maxPart];
  float STT, W, Yb;
  float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
    SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
    EC_Time[maxPart], EC_Path[maxPart],
    Charge[maxPart], Beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
    pz[maxPart], Theta[maxPart], Phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
    EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], EC_U[maxPart],EC_V[maxPart],EC_W[maxPart], 
    CC_Chi2[maxPart];
  t.SetBranchAddress("NRun"     ,&NRun   ); // Run number
  t.SetBranchAddress("gPart"    ,&gPart  ); // Number of particles observed (globally) in the event
  t.SetBranchAddress("CCPart"   ,&CCPart ); // Number of particles observed in the Cherenkovs
  t.SetBranchAddress("DCPart"   ,&DCPart ); // Number of particles observed in the Drift Chambers
  t.SetBranchAddress("ECPart"   ,&ECPart ); // Number of particles observed in the ECal
  t.SetBranchAddress("SCPart"   ,&SCPart ); // Number of particles observed in the Scintillation Counters (ToFs)
  t.SetBranchAddress("Stat"     ,Stat    ); // Global status for each particle candidate
  t.SetBranchAddress("StatDC"   ,StatDC  ); // Drift Chamber status for each particle candidate
  t.SetBranchAddress("StatCC"   ,StatCC  ); // Cherenkov status for each particle
  t.SetBranchAddress("StatEC"   ,StatEC  ); // ECal status for each particle
  t.SetBranchAddress("StatSC"   ,StatSC  ); // Scintillation counter status for each particle  
  t.SetBranchAddress("particle" ,id_guess); // Guess at the particle ID made by the recon software (maybe not reliable)
  t.SetBranchAddress("EC_in"    ,EC_in   ); // Inner layer of ECal for each particle
  t.SetBranchAddress("EC_out"   ,EC_out  ); // Outer layer of ECal for each particle
  t.SetBranchAddress("EC_tot"   ,EC_tot  ); // Total energy deposit in the ECal for each particle
  t.SetBranchAddress("Nphe"     ,Nphe    ); // Number of photo-electrons per hit in the Cherenkov detectors
  t.SetBranchAddress("SC_Time"  ,SC_Time ); // Time in the scintillators per particle
  t.SetBranchAddress("SC_Path"  ,SC_Path ); // Path Length per particle
  t.SetBranchAddress("SC_Pad"   ,SC_Pad );
  t.SetBranchAddress("CC_Time"  ,CC_Time ); // Time in the cherenkov per particle
  t.SetBranchAddress("CC_Path"  ,CC_Path ); // Path Length per particle
  t.SetBranchAddress("EC_Time"  ,EC_Time ); // Time in the EC per particle
  t.SetBranchAddress("EC_Path"  ,EC_Path ); // Path Length per particle
  t.SetBranchAddress("Charge"   ,Charge  ); // Charge per particle
  t.SetBranchAddress("Beta"     ,Beta    ); // Beta per particle
  t.SetBranchAddress("Mass"     ,mass    ); // Mass per particle
  t.SetBranchAddress("Momentum" ,mom     ); // Momentum magnitude per particle
  t.SetBranchAddress("Momentumx",px      ); // Momentum x component per particle
  t.SetBranchAddress("Momentumy",py      ); // Momentum y component per particle
  t.SetBranchAddress("Momentumz",pz      ); // Momentum z component per particle
  t.SetBranchAddress("Theta"    ,Theta   ); // Theta per particle
  t.SetBranchAddress("Phi"      ,Phi     ); // Phi per particle
  t.SetBranchAddress("TargetZ"  ,targetZ ); // Target Z per particle
  t.SetBranchAddress("Thetapq"  ,theta_pq); // Angle wrt to q vector per particle
  t.SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
  t.SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
  t.SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
  t.SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter
  t.SetBranchAddress("EC_U"     ,EC_U    ); // u positions of hit in the calorimeter
  t.SetBranchAddress("EC_V"     ,EC_V    ); // v positions of hit in the calorimeter
  t.SetBranchAddress("EC_W"     ,EC_W    ); // w positions of hit in the calorimeter
  t.SetBranchAddress("CC_Chi2"  ,CC_Chi2 ); // angle between CC hit and nearest SC hit (in rad)
  //t.SetBranchAddress("W"        ,&W      ); // Hadronic mass
  //t.SetBranchAddress("Yb"       ,&Yb     ); // Y-scaling variable

  /*
     the kinematic variables: Q2, W, Xb, Nu, Yb are calculated incorrectly in the particle_data.root file
     (since they assume the eg2c beam energy as hard-coded in the Tidentificator library).
     So, need to calculate them by hand.
  */

  // --------------------------------------------------------------------------------------------------
  // Open up the output file
  TFile * outfile = new TFile(argv[1],"RECREATE");
  TDirectory * dir_diag = outfile->mkdir("diagnostics");
  dir_diag->cd();

  // ---------------------------------------
  // Diagnostic electron histograms
  TDirectory * dir_e = dir_diag->mkdir("e");
  TH1D * h1_Xb0        = new TH1D("h1_Xb0"        ,"x_{B} without e- p correction"                     ,200,  0.8, 2.0);
  TH1D * h1_Xb1        = new TH1D("h1_Xb1"        ,"x_{B} with    e- p correction"                     ,200,  0.8, 2.0);
  
  TH1D * h1_e_Nphe0    = new TH1D("h1_e_Nphe0"    ,"e- before cuts;# photo-electrons in CC;Counts"     ,200,   0.,200.);
  TH1D * h1_e_Nphe1    = new TH1D("h1_e_Nphe1"    ,"e- passing PID;# photo-electrons in CC;Counts"     ,200,   0.,200.);
  TH1D * h1_e_Nphe2    = new TH1D("h1_e_Nphe2"    ,"e- passing PID+fid;# photo-electrons in CC;Counts" ,200,   0.,200.);
  
  TH1D * h1_e_EC_in0   = new TH1D("h1_e_EC_in0"   ,"e- before cuts;E_in [GeV];Counts"                  ,300,   0.,  1.);
  TH1D * h1_e_EC_in1   = new TH1D("h1_e_EC_in1"   ,"e- passing PID;E_in [GeV];Counts"                  ,300,   0.,  1.);
  TH1D * h1_e_EC_in2   = new TH1D("h1_e_EC_in2"   ,"e- passing PID+fid;E_in [GeV];Counts"              ,300,   0.,  1.);
  
  TH1D * h1_e_EC_out0  = new TH1D("h1_e_EC_out0"  ,"e- before cuts;E_out [GeV];Counts"                 ,300,   0.,  1.);
  TH1D * h1_e_EC_out1  = new TH1D("h1_e_EC_out1"  ,"e- passing PID;E_out [GeV];Counts"                 ,300,   0.,  1.);
  TH1D * h1_e_EC_out2  = new TH1D("h1_e_EC_out2"  ,"e- passing PID+fid;E_out [GeV];Counts"             ,300,   0.,  1.);
  
  TH1D * h1_e_EC_tot0   = new TH1D("h1_e_EC_tot0"  ,"e- before cuts;E_tot [GeV];Counts"                ,300,   0.,  1.);
  TH1D * h1_e_EC_tot1   = new TH1D("h1_e_EC_tot1"  ,"e- passing PID;E_tot [GeV];Counts"                ,300,   0.,  1.);
  TH1D * h1_e_EC_tot2   = new TH1D("h1_e_EC_tot2"  ,"e- passing PID+fid;E_tot [GeV];Counts"            ,300,   0.,  1.);
  
  TH2D * h2_e_phiTheta0= new TH2D("h2_e_phiTheta0","e- before  cuts;#phi [deg];#theta [deg];Counts"      ,300,-100.,380.,300,10.,50.);
  TH2D * h2_e_phiTheta1= new TH2D("h2_e_phiTheta1","e- passing PID;#phi [deg];#theta [deg];Counts"       ,300,-100.,380.,300,10.,50.);
  TH2D * h2_e_phiTheta2= new TH2D("h2_e_phiTheta2","e- passing PID+fid;#phi [deg];#theta [deg];Counts"   ,300,-100.,380.,300,10.,50.);
  
  TH2D * h2_e_Ein_Eout0= new TH2D("h2_e_Ein_Eout0","e- before cuts;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
  TH2D * h2_e_Ein_Eout1= new TH2D("h2_e_Ein_Eout1","e- passing PID;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
  TH2D * h2_e_Ein_Eout2= new TH2D("h2_e_Ein_Eout2","e- passing PID+fid;E_in/p;E_out/p;Counts"          ,300,   0., 0.5,300, 0.,0.5);
  
  TH2D * h2_e_Ein_Eout_0=new TH2D("h2_e_Ein_Eout_0","e- before cuts;E_in;E_out;Counts"                 ,300,   0., 1.5,300, 0.,1.5);
  TH2D * h2_e_Ein_Eout_1=new TH2D("h2_e_Ein_Eout_1","e- passing PID;E_in;E_out;Counts"                 ,300,   0., 1.5,300, 0.,1.5);
  TH2D * h2_e_Ein_Eout_2=new TH2D("h2_e_Ein_Eout_2","e- passing PID+fid;E_in;E_out;Counts"             ,300,   0., 1.5,300, 0.,1.5);
  
  TH2D * h2_e_xyEC_hit0= new TH2D("h2_e_xyEC_hit0","e- before cuts;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
  TH2D * h2_e_xyEC_hit1= new TH2D("h2_e_xyEC_hit1","e- passing PID;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
  TH2D * h2_e_xyEC_hit2= new TH2D("h2_e_xyEC_hit2","e- passing PID+fid;ECx [cm];ECy [cm];Counts"       ,300,-400.,400.,300,-400.,400.);
  
  TH2D * h2_e_p_Etot0  = new TH2D("h2_e_p_Etot0"  ,"e- before cuts;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
  TH2D * h2_e_p_Etot1  = new TH2D("h2_e_p_Etot1"  ,"e- passing PID;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
  TH2D * h2_e_p_Etot2  = new TH2D("h2_e_p_Etot2"  ,"e- passing PID+fid;p [GeV];E_tot/p;Counts"         ,300,   0.,  5.,300, 0.,0.7);
  
  TH2D * h2_e_p_E0     = new TH2D("h2_e_p_E0"     ,"e- before cuts;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
  TH2D * h2_e_p_E1     = new TH2D("h2_e_p_E1"     ,"e- passing PID;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
  TH2D * h2_e_p_E2     = new TH2D("h2_e_p_E2"     ,"e- before PIF+fid;p [GeV];E;Counts"                ,300,   0.,  5.,300, 0.,2.);
  
  TH2D * h2_e_thetaMom0= new TH2D("e_thetaMom0"     ,"e- before cuts;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
  TH2D * h2_e_thetaMom1= new TH2D("e_thetaMom1"     ,"e- passing PID;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
  TH2D * h2_e_thetaMom2= new TH2D("e_thetaMom2"     ,"e- passing PID+fid;#theta [deg];Mom [GeV];Counts"   ,300,  10., 50.,300, 0., 6.);
  
  // ---
  TH1D * h1_e_vz_sec10 = new TH1D("h1_e_vz_sec10" ,"e- passing cuts, before vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec1  = new TH1D("h1_e_vz_sec1"  ,"e- passing cuts,  after vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec20 = new TH1D("h1_e_vz_sec20" ,"e- passing cuts, before vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec2  = new TH1D("h1_e_vz_sec2"  ,"e- passing cuts,  after vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec30 = new TH1D("h1_e_vz_sec30" ,"e- passing cuts, before vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec3  = new TH1D("h1_e_vz_sec3"  ,"e- passing cuts,  after vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec40 = new TH1D("h1_e_vz_sec40" ,"e- passing cuts, before vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec4  = new TH1D("h1_e_vz_sec4"  ,"e- passing cuts,  after vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec53 = new TH1D("h1_e_vz_sec53" ,"e- passing cuts, before vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec5  = new TH1D("h1_e_vz_sec5"  ,"e- passing cuts,  after vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec60 = new TH1D("h1_e_vz_sec60" ,"e- passing cuts, before vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz_sec6  = new TH1D("h1_e_vz_sec6"  ,"e- passing cuts,  after vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
  TH1D * h1_e_vz0      = new TH1D("h1_e_vz0"      ,"e- passing cuts, before vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
  TH1D * h1_e_vz       = new TH1D("h1_e_vz"       ,"e- passing cuts,  after vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
  TH2D * h2_e_phiVz0   = new TH2D("h2_e_phiVz0"   ,"e- passing cuts, before vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
  TH2D * h2_e_phiVz    = new TH2D("h2_e_phiVz"    ,"e- passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
  TH2D * h2_e_thetaVz0 = new TH2D("h2_e_thetaVz0" ,"e- passing cuts, before vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
  TH2D * h2_e_thetaVz  = new TH2D("h2_e_thetaVz"  ,"e- passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
  // ---
  
  TH1D * h1_e_momCor    = new TH1D("h1_e_momCor"    ,"e- passing fid. cuts;p corrected - p[GeV];Counts"          ,300, -.1,.04 );
  TH1D * h1_e_momCor1   = new TH1D("h1_e_momCor1"   ,"e- passing fid. cuts;p corrected/p;Counts"                 ,300,0.97,1.01);
  TH2D * h2_e_momMomCor = new TH2D("h2_e_momMomCor" ,"e- passing fid. cuts;p [GeV];p corrected - p[GeV];Counts"  ,300,  0.,  6.,300,-.1 ,.04 );
  TH2D * h2_e_momMomCor1= new TH2D("h2_e_momMomCor1","e- passing fid. cuts;p [GeV];p corrected/p;Counts"         ,300,  0.,  6.,300,0.97,1.01);

  TH1D * h1_e_momCor_sec1   = new TH1D("h1_e_momCor_sec1"   ,"e- passing PID+fid sec1;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec1= new TH2D("h2_e_momMomCor_sec1","e- passing PID+fid sec1;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  TH1D * h1_e_momCor_sec2   = new TH1D("h1_e_momCor_sec2"   ,"e- passing PID+fid sec2;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec2= new TH2D("h2_e_momMomCor_sec2","e- passing PID+fid sec2;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  TH1D * h1_e_momCor_sec3   = new TH1D("h1_e_momCor_sec3"   ,"e- passing PID+fid sec3;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec3= new TH2D("h2_e_momMomCor_sec3","e- passing PID+fid sec3;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  TH1D * h1_e_momCor_sec4   = new TH1D("h1_e_momCor_sec4"   ,"e- passing PID+fid sec4;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec4= new TH2D("h2_e_momMomCor_sec4","e- passing PID+fid sec4;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  TH1D * h1_e_momCor_sec5   = new TH1D("h1_e_momCor_sec5"   ,"e- passing PID+fid sec5;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec5= new TH2D("h2_e_momMomCor_sec5","e- passing PID+fid sec5;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  TH1D * h1_e_momCor_sec6   = new TH1D("h1_e_momCor_sec6"   ,"e- passing PID+fid sec6;p corrected/p;Counts"        ,300,0.97,1.01);
  TH2D * h2_e_momMomCor_sec6= new TH2D("h2_e_momMomCor_sec6","e- passing PID+fid sec6;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
  
  TH2D * h2_e_vzVzCor   = new TH2D("h2_e_vzVzCor"   ,"e- passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300,-20., 20.,300,-1. , 1. );
  
  // Diagnostic electron momentum correction
  TH2D * h2_e_thetaMom3= new TH2D("e_thetaMom3","e- passing PID+fid (z axis: p_{corr}-p);#theta [deg];Mom [GeV];p correction" ,300,  10., 50.,300, 0., 6.);

  // ---------------------------------------
  // Diagnostic positive particle histograms
  TDirectory * dir_pos = dir_diag->mkdir("pos");
  TH1D * h1_p_mass     = new TH1D("h1_pos_mass"   ,"+  passing fid. cuts;mass [GeV];Counts"          ,300,  0.,3.5 );
  TH2D * h2_p_pMass    = new TH2D("h2_pos_pMass"  ,"+  passing fid. cuts;p [GeV];mass [GeV];Counts"  ,300,  0., 4.,300, 0.,3.5);
  TH2D * h2_pos_pBeta  = new TH2D("h2_pos_pBeta"  ,"+  passing fid. cuts;p [GeV];#beta;Counts"       ,300,  0., 4.,300, 0.,1.3);
  
  // ---------------------------------------
  // Diagnostic proton histograms
  TDirectory * dir_p = dir_diag->mkdir("p");
  TH2D * h2_p_phiTheta0= new TH2D("h2_p_phiTheta0","p before cuts;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
  TH2D * h2_p_phiTheta1= new TH2D("h2_p_phiTheta1","p passing PID;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
  TH2D * h2_p_phiTheta2= new TH2D("h2_p_phiTheta2","p passing fid+PID;#phi [deg];#theta [deg];Counts",300,-100.,380.,300,0.,55.);
  
  TH2D * h2_p_deltaTmom0=new TH2D("h2_p_deltaTmom0","p before cuts;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 6.);
  TH2D * h2_p_deltaTmom1=new TH2D("h2_p_deltaTmom1","p passing PID;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 6.);
  TH2D * h2_p_deltaTmom2=new TH2D("h2_p_deltaTmom2","p passing fid+PID;#Delta t [ns];p [GeV];Counts" ,300,  -7.,  7.,300, 0., 6.);
  
  TH2D * h2_p_p_momCor0= new TH2D("h2_p_p_momCor0","p passing fid. cuts;p [GeV];p - p_corr [GeV];Counts"     ,100, 0., 2.5,100,-.06,.01);
  TH2D * h2_p_p_momCor1= new TH2D("h2_p_p_momCor1","p passing fid. cuts;p [GeV];p_corr/p;Counts"             ,100, 0., 2.5,100,0.95,1.2);
  TH2D * h2_p_th_pCor0 = new TH2D("h2_p_th_pCor0" ,"p passing fid. cuts;#theta [deg];p - p_corr [GeV];Counts",100, 0., 150,100,-.06,.01);
  TH2D * h2_p_th_pCor1 = new TH2D("h2_p_th_pCor1" ,"p passing fid. cuts;#theta [deg];p_corr/p;Counts"        ,100, 0., 150,100,0.95,1.2);
  TH2D * h2_p_th_p_cor0= new TH2D("h2_p_th_p_cor0","p passing fid. cuts;#theta [deg];p[GeV];p - p_corr [GeV]",100, 0., 150,100,0.  ,2.5);
  TH2D * h2_p_th_p_cor1= new TH2D("h2_p_th_p_cor1","p passing fid. cuts;#theta [deg];p[GeV];p_corr/p"        ,100, 0., 150,100,0.  ,2.5);
  
  TH2D * h2_p_vzVzCor  = new TH2D("h2_p_vzVzCor"  ,"p passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300, -20., 20.,300,-1., 1 );
  TH2D * h2_p_phiVz0   = new TH2D("h2_p_phiVz0"   ,"p passing cuts, before vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
  TH2D * h2_p_phiVz    = new TH2D("h2_p_phiVz"    ,"p passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
  TH2D * h2_p_thetaVz0 = new TH2D("h2_p_thetaVz0" ,"p passing cuts, before vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
  TH2D * h2_p_thetaVz  = new TH2D("h2_p_thetaVz"  ,"p passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
  TH2D * h2_p_pBeta    = new TH2D("h2_p_pBeta"    ,"p passing fid. cuts;p [GeV];#beta;Counts"                 ,300,   0.,  4.,300, 0.,1.3);
  
  // ---------------------------------------
  // Diagnostic pi+ histograms
  TDirectory * dir_pip = dir_diag->mkdir("pi+");
  TH2D * h2_pip_pBeta     = new TH2D("h2_pip_pBeta"     ,"#pi+ passing fid. cuts;p [GeV];#beta;Counts"      ,300,   0.,  4.,300, 0.,1.3);
  TH2D * h2_pip_deltaTmom0= new TH2D("h2_pip_deltaTmom0","#pi+ before cuts;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
  TH2D * h2_pip_deltaTmom1= new TH2D("h2_pip_deltaTmom1","#pi+ passing fid;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
  TH2D * h2_pip_deltaTmom2= new TH2D("h2_pip_deltaTmom2","#pi+ passing fid+PID;#Delta t [ns];p [GeV];Counts",300,  -7.,  7.,300, 0.,5.0);
  
  // ---------------------------------------
  // Diagnostic neutral histograms
  TDirectory * dir_neut = dir_diag->mkdir("neut");  
  TH2D * h2_neu_pBeta  = new TH2D("h2_neu_pBeta"  ,"neutral passing fid. cuts;p [GeV];#beta;Counts"   ,300,   0.,  4.,300, 0.,1.3);
  
  // ---------------------------------------
  // Diagnostic neutron histograms
  TDirectory * dir_n = dir_diag->mkdir("n");  
  TH2D * h2_n_phiTheta0= new TH2D("h2_n_phiTheta0","n before cuts;#phi [deg];#theta [deg];Counts"     ,300,-100.,380.,300,0.,55.);
  TH2D * h2_n_phiTheta1= new TH2D("h2_n_phiTheta1","n passing fid;#phi [deg];#theta [deg];Counts"     ,300,-100.,380.,300,0.,55.);
  TH2D * h2_n_phiTheta2= new TH2D("h2_n_phiTheta2","n passing fid+PID;#phi [deg];#theta [deg];Counts" ,300,-100.,380.,300,0.,55.);
  
  TH2D * h2_n_pBeta    = new TH2D("h2_n_pBeta"    ,"n passing fid. cuts;p [GeV];#beta;Counts"         ,300,   0.,  4.,300, 0.,1.3);
  
  TH1D * h1_u_0        = new TH1D("h1_u_0"        ,"n before xyz cut;EC_{u} [cm];Counts"              ,100,   0., 500.);
  TH1D * h1_v_0        = new TH1D("h1_v_0"        ,"n before xyz cut;EC_{v} [cm];Counts"              ,100,   0., 500.);
  TH1D * h1_w_0        = new TH1D("h1_w_0"        ,"n before xyz cut;EC_{w} [cm];Counts"              ,100,   0., 500.);
  TH1D * h1_u_1        = new TH1D("h1_u_1"        ,"n after xyz cut;EC_{u} [cm];Counts"               ,100,   0., 500.);
  TH1D * h1_v_1        = new TH1D("h1_v_1"        ,"n after xyz cut;EC_{v} [cm];Counts"               ,100,   0., 500.);
  TH1D * h1_w_1        = new TH1D("h1_w_1"        ,"n after xyz cut;EC_{w} [cm];Counts"               ,100,   0., 500.);
  
  TH1D * h1_x_0        = new TH1D("h1_x_0"        ,"n before xyz cut;EC_{x} [cm];Counts"              ,100, -500., 500.);
  TH1D * h1_y_0        = new TH1D("h1_y_0"        ,"n before xyz cut;EC_{y} [cm];Counts"              ,100, -500., 500.);
  TH1D * h1_z_0        = new TH1D("h1_z_0"        ,"n before xyz cut;EC_{z} [cm];Counts"              ,100,  300., 600.);
  TH1D * h1_x_1        = new TH1D("h1_x_1"        ,"n after xyz cut;EC_{x} [cm];Counts"               ,100, -500., 500.);
  TH1D * h1_y_1        = new TH1D("h1_y_1"        ,"n after xyz cut;EC_{y} [cm];Counts"               ,100, -500., 500.);
  TH1D * h1_z_1        = new TH1D("h1_z_1"        ,"n after xyz cut;EC_{z} [cm];Counts"               ,100,  300., 600.);
  
  TH2D * h2_n_ECxy_0   = new TH2D("h2_n_ECxy_0"   ,"n before xyz cut;EC_{x} [cm];EC_{y} [cm];Counts"  ,300, -500., 500., 300, -500., 500.);
  TH2D * h2_n_ECxy_1   = new TH2D("h2_n_ECxy_1"   ,"n after xyz cut;EC_{x} [cm];EC_{y} [cm];Counts"   ,300, -500., 500., 300, -500., 500.);

  // ---------------------------------------
  // Diagnostic pi- histograms
  TDirectory * dir_pim = dir_diag->mkdir("pi-");
  TH2D * h2_pim_pBeta     = new TH2D("h2_pim_pBeta"     ,"#pi- passing fid. cuts;p [GeV];#beta;Counts"      ,300,   0.,  4.,300, 0.,1.3);
  TH2D * h2_pim_deltaTmom0= new TH2D("h2_pim_deltaTmom0","#pi- before cuts;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
  TH2D * h2_pim_deltaTmom1= new TH2D("h2_pim_deltaTmom1","#pi- passing fid;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
  TH2D * h2_pim_deltaTmom2= new TH2D("h2_pim_deltaTmom2","#pi- passing fid+PID;#Delta t [ns];p [GeV];Counts",300,  -7.,  7.,300, 0.,5.0);
  
  // ---------------------------------------
  // Setting up output tree and branches
  TTree * outtree = new TTree("T","Skimmed tree");
  double e_vz, e_vz_corrected, e_mom[3], e_phi_mod;
  double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
  double e_t0,beta_assuming_proton,p_t0,delta_t,beta_assuming_pip,pip_t0,pip_delta_t;
  double corr_px, corr_py, corr_pz, n_px, n_py, n_pz, n_p, EC_Path_corr, Beta_corr;
  
  int nRun, nParticles;
  int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
  int Part_type    [maxPart];
  double Nu, Q2, Xb, Nu_unc, Q2_unc, Xb_unc, t0;
  double vtx_z_unc [maxPart], vtx_z_cor[maxPart], Mass[maxPart];
  double mom_x     [maxPart], mom_y    [maxPart], mom_z  [maxPart];
  double theta[maxPart], phi[maxPart];
  double e_deltat  [maxPart];
  int    stat_sc   [maxPart], stat_ec  [maxPart], stat_dc[maxPart], sc_pad[maxPart];
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

  outtree->Branch("nRun"      , &nRun      , "nRun/I"                  );
  outtree->Branch("nParticles", &nParticles, "nParticles/I"            );
  outtree->Branch("nProtons"  , &nProtons  , "nProtons/I"              );
  outtree->Branch("nNeutrons" , &nNeutrons , "nNeutrons/I"             );
  outtree->Branch("nPiplus"   , &nPiplus   , "nPiplus/I"               );
  outtree->Branch("nPiminus"  , &nPiminus  , "nPiminus/I"              );
  outtree->Branch("t0"        , &t0        , "t0/D"                    );
  outtree->Branch("Nu"        , &Nu        , "Nu/D"                    );
  outtree->Branch("Q2"        , &Q2        , "Q2/D"                    );
  outtree->Branch("Xb"        , &Xb        , "Xb/D"                    );
  outtree->Branch("charge"    ,  charge    , "charge[nParticles]/D"    );
  outtree->Branch("beta"      ,  beta      , "beta[nParticles]/D"      );
  outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" );
  outtree->Branch("vtx_z_unc" ,  vtx_z_unc , "vtx_z_unc[nParticles]/D" );
  outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" );
  outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     );
  outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     );
  outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     );
  outtree->Branch("phi"     ,  phi     , "phi[nParticles]/D"     );
  outtree->Branch("theta"     ,  theta     , "theta[nParticles]/D"     );
  outtree->Branch("e_deltat"  ,  e_deltat  , "e_deltat[nParticles]/D"  );
  outtree->Branch("stat_sc"   ,  stat_sc   , "stat_sc[nParticles]/I"   );
  outtree->Branch("stat_dc"   ,  stat_dc   , "stat_dc[nParticles]/I"   );
  outtree->Branch("stat_ec"   ,  stat_ec   , "stat_ec[nParticles]/I"   );
  outtree->Branch("sc_time"   ,  sc_time   , "sc_time[nParticles]/D"   );
  outtree->Branch("sc_path"   ,  sc_path   , "sc_path[nParticles]/D"   );
  outtree->Branch("sc_pad"    ,  sc_pad    , "sc_pad[nParticles]/I"    );
  outtree->Branch("ec_time"   ,  ec_time   , "ec_time[nParticles]/D"   );
  outtree->Branch("ec_path"   ,  ec_path   , "ec_path[nParticles]/D"   );
  outtree->Branch("ec_in"     ,  ec_in     , "ec_in[nParticles]/D"     );
  outtree->Branch("ec_out"    ,  ec_out    , "ec_out[nParticles]/D"    );
  outtree->Branch("ec_tot"    ,  ec_tot    , "ec_tot[nParticles]/D"    );
  outtree->Branch("ec_x"      ,  ec_x      , "ec_x[nParticles]/D"      );
  outtree->Branch("ec_y"      ,  ec_y      , "ec_y[nParticles]/D"      );
  outtree->Branch("ec_z"      ,  ec_z      , "ec_z[nParticles]/D"      );
  outtree->Branch("ec_u"      ,  ec_u      , "ec_u[nParticles]/D"      );
  outtree->Branch("ec_v"      ,  ec_v      , "ec_v[nParticles]/D"      );
  outtree->Branch("ec_w"      ,  ec_w      , "ec_w[nParticles]/D"      );
  outtree->Branch("Mass"      ,  Mass      , "Mass[nParticles]/D"      );
  // --------------------------------------------------------------------------------------------------
  // Obtaining run number and other important parameters
  t.GetEvent(0);

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
  } while(tab_run != NRun);
  
  tab_run=NRun;
  Run_dependent run_dependent_corrections(NRun);       // Create an instance of the Run_dependent Class
  tab_E1    = run_dependent_corrections.get_E1();
  tab_torus = run_dependent_corrections.get_torus();
  tab_mini  = run_dependent_corrections.get_mini();
  tab_targ  = run_dependent_corrections.get_targ();

  cout << "Run    = " << tab_run   << endl;
  cout << "Ebeam  = " << tab_E1    << endl;
  cout << "Torus  = " << tab_torus << endl;
  cout << "Mini   = " << tab_mini  << endl;
  cout << "Target = " << tab_targ  << endl;  
  Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, true);  // Create an instance of the Fiducial Class

  // Values for some cuts
  const double EC_in_cut = fid_params.EC_in_cut();
  const double el_EC_cut = fid_params.el_EC_cut();
  
  // --------------------------------------------------------------------------------------------------
  // Loop over events
  for (int event=0; event < nEvents ; event++)

    {
      if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}

      t.GetEvent(event);
      nRun = NRun;

      if (gPart <= 0) continue; // Ignore events that have no particle candidates

      // Set all the counters to zero for the start of the event
      nParticles = 0;
      nProtons   = 0;
      nNeutrons  = 0;
      nPiplus    = 0;
      nPiminus   = 0;
      nPi0       = 0;

      // --------------------------------------------------------------------------------------------------
      // Sector index for electrons
      int e_sect = (int)(Phi[0]+30)/60;
      if (e_sect>5) e_sect = 5;
      if (e_sect<0) e_sect = 0;

      // --------------------------------------------------------------------------------------------------
      double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC
      TVector3 T3_e_mom(px[0],py[0],pz[0]); // Electron momentum expressed in a TVector3
      e_vz_corrected = targetZ[0]+fid_params.vz_corr(T3_e_mom);
      TVector3 e_ec_xyz(EC_X[0],EC_Y[0],EC_Z[0]);

      // ---------------------------------------------------------------------------------------
      // Electron general cuts
      if (!(			(StatEC[0] > 0) && // EC status is good for the electron candidate
                  (StatDC[0] > 0) && // DC status is good for the electron candidate
                  (StatCC[0] > 0) && // CC status is good for the electron candidate
                  (StatSC[0] > 0) && // SC status is good for the electron candidate
                  (Charge[0] < 0)    // Electron candidate curvature direction is negative
                  ))
        {continue;}
      // ---------------------------------------------------------------------------------------

      h1_e_Nphe0     -> Fill( Nphe        [0]  );
      h1_e_EC_in0    -> Fill( EC_in       [0]  );
      h1_e_EC_out0   -> Fill( EC_out      [0]  );
      h1_e_EC_tot0   -> Fill( EC_tot      [0]  );
      h2_e_phiTheta0 -> Fill( Phi         [0], Theta        [0]);
      h2_e_Ein_Eout0 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
      h2_e_Ein_Eout_0-> Fill( EC_in[0]       , EC_out[0]       );
      h2_e_p_Etot0   -> Fill( mom         [0], EC_tot[0]/mom[0]);
      h2_e_xyEC_hit0 -> Fill( EC_X        [0], EC_Y         [0]);
      h2_e_p_E0      -> Fill( mom         [0], el_cand_EC      );
      h2_e_thetaMom0 -> Fill( Theta       [0], mom          [0]);

      // ---------------------------------------------------------------------------------------
      //Electron particle Identification
      if (!(                  (EC_in [0] > EC_in_cut) &&      // Electron candidate has enough energy deposit in inner layer of EC    
                              (el_cand_EC > el_EC_cut) &&     // Enough total energy in the EC
                              (fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) // Electron PID (E/p)
                              ))
        {continue;}
      
      // ---------------------------------------
      // Additional cut for 2GeV data:
      double el_sccc_dt = SC_Time[0] - CC_Time[0] - (SC_Path[0] - CC_Path[0])/(c_m_s*ns_to_s*100.);
      
      if(			(tab_E1==2261)&&(
                               CC_Chi2[0]>=0.1 ||
                               el_sccc_dt < sc_cc_dt_cut_sect[e_sect] ||
                               sqrt(mom[0]*mom[0]+me*me)>tab_E1/1000.
                               ))
        {continue;}

      // ---------------------------------------------------------------------------------------
      // If the event made it here, the electron candidate passed all PID cuts
      h1_e_Nphe1     -> Fill( Nphe        [0]  );
      h1_e_EC_in1    -> Fill( EC_in       [0]  );
      h1_e_EC_out1   -> Fill( EC_out      [0]  );
      h1_e_EC_tot1   -> Fill( EC_tot      [0]  );
      h2_e_phiTheta1 -> Fill( Phi         [0], Theta        [0]);
      h2_e_Ein_Eout1 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
      h2_e_Ein_Eout_1-> Fill( EC_in[0]       , EC_out[0]       );
      h2_e_p_Etot1   -> Fill( mom         [0], EC_tot[0]/mom[0]);
      h2_e_xyEC_hit1 -> Fill( EC_X        [0], EC_Y         [0]);
      h2_e_p_E1      -> Fill( mom         [0], el_cand_EC      );
      h2_e_thetaMom1 -> Fill( Theta       [0], mom          [0]);

      // ---------------------------------------------------------------------------------------
      //cout << phi[0] << endl;
      if (do_fid_e)
        {
          // Electron Fiducial cuts
          if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
          if (!fid_params.CutUVW_e(e_ec_xyz)       ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);
        }
      //cout << phi[0] << endl;
      // ---------------------------------------------------------------------------------------

      // If electron passes all cuts, then apply momentum corrections
      //TVector3 T3_e_mom_cor = T3_e_mom; // This skips momentum correction
      TVector3 T3_e_mom_cor = fid_params.eMomentumCorrection(T3_e_mom);
      
      // If we get to here, then the electron passed fiducial cuts
      t0 = STT;

      // Without electron momentum correction
      Nu_unc = tab_E1/1000. - T3_e_mom.Mag(); //Energy transfer
      Q2_unc = 4.*tab_E1/1000.*T3_e_mom.Mag()*sin(T3_e_mom.Theta()/2.)*sin(T3_e_mom.Theta()/2.);      //4-momentum transfer^2
      Xb_unc = Q2_unc / (2*mP*Nu_unc);    //Bjorken scaling variable

      // With electron momentum correction (these are saved in the tree)
      Nu = tab_E1/1000. - T3_e_mom_cor.Mag();	//Energy transfer
      Q2 = 4.*tab_E1/1000.*T3_e_mom_cor.Mag()*sin(T3_e_mom_cor.Theta()/2.)*sin(T3_e_mom_cor.Theta()/2.);	//4-momentum transfer^2
      Xb = Q2 / (2*mP*Nu);	//Bjorken scaling variable

      nParticles++;
      Part_type[0] = -11;
      e_deltat [0] =   0;
      mom_x    [0] = T3_e_mom_cor.X();
      mom_y    [0] = T3_e_mom_cor.Y();
      mom_z    [0] = T3_e_mom_cor.Z();
      vtx_z_unc[0] = targetZ[0];
      vtx_z_cor[0] = e_vz_corrected;
      stat_sc  [0] = StatSC [0];
      stat_dc  [0] = StatDC [0];
      stat_ec  [0] = StatEC [0];
      sc_time  [0] = SC_Time[0];
      sc_path  [0] = SC_Path[0];
      sc_pad  [0] = SC_Pad[0];
      ec_time  [0] = EC_Time[0];
      ec_path  [0] = EC_Path[0];
      ec_in    [0] = EC_in  [0];
      ec_out   [0] = EC_out [0];
      ec_tot   [0] = EC_tot [0];
      ec_x     [0] = EC_X   [0];
      ec_y     [0] = EC_Y   [0];
      ec_z     [0] = EC_Z   [0];
      ec_u     [0] = EC_U   [0];
      ec_v     [0] = EC_V   [0];
      ec_w     [0] = EC_W   [0];
      Mass     [0] = mass   [0];
      charge   [0] = Charge [0];
      beta     [0] = Beta   [0];
      theta[0] = Theta[0];
      phi[0] = Phi[0];

      // Fill some diagnostic histograms
      h1_Xb0         -> Fill(Xb_unc   );
      h1_Xb1         -> Fill(Xb       );
      h1_e_Nphe2     -> Fill(Nphe  [0]);
      h1_e_EC_in2    -> Fill(EC_in [0]);
      h1_e_EC_out2   -> Fill(EC_out[0]);
      h1_e_EC_tot2   -> Fill(EC_tot[0]);
      h2_e_xyEC_hit2 -> Fill(EC_X[0],EC_Y[0]);
      h2_e_thetaMom2 -> Fill(Theta[0],mom[0]);
      h1_e_momCor    -> Fill(T3_e_mom_cor.Mag()-T3_e_mom.Mag());
      h1_e_momCor1   -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      h2_e_momMomCor -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()-T3_e_mom.Mag());
      h2_e_momMomCor1-> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      h2_e_vzVzCor   -> Fill(targetZ[0],e_vz_corrected-targetZ[0]);
      h2_e_Ein_Eout2 -> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
      h2_e_Ein_Eout_2-> Fill(EC_in  [0]     , EC_out[0]      );
      h2_e_p_Etot2   -> Fill(mom    [0],EC_tot[0]/mom[0]);
      h2_e_p_E2      -> Fill(mom    [0],el_cand_EC      );
      h2_e_phiTheta2 -> Fill(Phi    [0],Theta  [0]      );
      h2_e_phiVz0    -> Fill(Phi    [0],targetZ[0]      );
      h2_e_phiVz     -> Fill(Phi    [0],e_vz_corrected  );
      h2_e_thetaVz0  -> Fill(Theta  [0],targetZ[0]      );
      h2_e_thetaVz   -> Fill(Theta  [0],e_vz_corrected  );
      h1_e_vz0       -> Fill(targetZ[0]                 );
      h1_e_vz        -> Fill(e_vz_corrected             );

      if      (e_sect==0) {
        h1_e_vz_sec10       -> Fill(targetZ[0]);
        h1_e_vz_sec1        -> Fill(e_vz_corrected);
        h1_e_momCor_sec1    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec1 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else if (e_sect==1) {
        h1_e_vz_sec20       -> Fill(targetZ[0]);
        h1_e_vz_sec2        -> Fill(e_vz_corrected);
        h1_e_momCor_sec2    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec2 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else if (e_sect==2) {
        h1_e_vz_sec30       -> Fill(targetZ[0]);
        h1_e_vz_sec3        -> Fill(e_vz_corrected);
        h1_e_momCor_sec3    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec3 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else if (e_sect==3) {
        h1_e_vz_sec40       -> Fill(targetZ[0]);
        h1_e_vz_sec4        -> Fill(e_vz_corrected);
        h1_e_momCor_sec4    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec4 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else if (e_sect==4) {
        h1_e_vz_sec53       -> Fill(targetZ[0]);
        h1_e_vz_sec5        -> Fill(e_vz_corrected);
        h1_e_momCor_sec5    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec5 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else if (e_sect==5) {
        h1_e_vz_sec60       -> Fill(targetZ[0]);
        h1_e_vz_sec6        -> Fill(e_vz_corrected);
        h1_e_momCor_sec6    -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
        h2_e_momMomCor_sec6 -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
      }
      else {cout << "Something is wrong with the definition of sectors" << endl;}
      h2_e_thetaMom3-> Fill(Theta[0],mom[0],T3_e_mom_cor.Mag()-T3_e_mom.Mag());
      // --------------------------------------------------------------------------------------------------
      // Loop over events looking for other particles

      for (int i=1 ; i<gPart ; i++)
        {
          Part_type[i] = 0; // Define the particle type as zero until a proper assignment can be made

          TVector3 T3_p_mom(px[i],py[i],pz[i]);
          TVector3 u1 = T3_p_mom.Unit();

          e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;

          beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
          p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
          delta_t = p_t0 - e_t0;

          beta_assuming_pip = mom[i]/sqrt(mom[i]*mom[i] + mpc*mpc);
          pip_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_pip * c_cm_ns);
          pip_delta_t = pip_t0 - e_t0;
          // ------------------------------------------------------------------------------------------
          // Test if positive particle
          if(             (StatSC[i] > 0) && 		// SC status is good for the positive candidate
                          (StatDC[i] > 0) &&              // DC status is good for the positive candidate
                          (Stat  [i] > 0) &&		// Global status is good for the positive candidate
                          (Charge[i] > 0) 		// Charge is positive
                          )
            {

              h2_p_phiTheta0    -> Fill(Phi[i]     ,Theta[i]);
              h2_p_deltaTmom0   -> Fill(delta_t    ,mom  [i]); 
              h2_pip_deltaTmom0 -> Fill(pip_delta_t,mom  [i]);
              
              // Positive particle vertex (_z) correction
              p_vz_corrected = targetZ[i]+fid_params.vz_corr(T3_p_mom);
              
              // --------------------------------------------------------------------
              // Look specifically for protons
              if( (max_p >= min_p) && // If we care about protons,
		  (id_guess[i] == 2212 ) &&  // ... and the particle is a proton candidate
		  (fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange))  // and the particle passes Proton PID cuts (delta T vs p)
                 ){
                h1_p_mass         -> Fill(mass[i]);
                h2_p_pMass        -> Fill(mom [i]    ,mass [i]);
                h2_pos_pBeta      -> Fill(mom [i]    ,Beta [i]);
                h2_p_phiTheta1    -> Fill(Phi [i]    ,Theta[i]);
                h2_p_deltaTmom1   -> Fill(delta_t    ,mom  [i]);
                
                // passed PID, so check fiducials
                if (!do_fid_p || (do_fid_p && fid_params.pFiducialCut(T3_p_mom)))
                  {

                    //p_mom_corrected=mom[i];  // Ignore momentum corrections
                    if(run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected) != -1)
                      p_mom_corrected=run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected);
                    else	p_mom_corrected=mom[i];	

                    corr_px = p_mom_corrected*u1.X();
                    corr_py = p_mom_corrected*u1.Y();
                    corr_pz = p_mom_corrected*u1.Z();
		    
                    T3_p_mom.SetXYZ(corr_px,corr_py,corr_pz);
		    
                    Part_type[nParticles] = 2212;
                    e_deltat [nParticles] = delta_t;
                    mom_x    [nParticles] = T3_p_mom.X();
                    mom_y    [nParticles] = T3_p_mom.Y();
                    mom_z    [nParticles] = T3_p_mom.Z();
                    vtx_z_unc[nParticles] = targetZ  [i];	
                    vtx_z_cor[nParticles] = p_vz_corrected;
                    stat_sc  [nParticles] = StatSC [i];
                    stat_dc  [nParticles] = StatDC [i];
                    stat_ec  [nParticles] = StatEC [i];
                    sc_time  [nParticles] = SC_Time[i];
                    sc_path  [nParticles] = SC_Path[i];
                    sc_pad  [nParticles] = SC_Pad[i];
                    ec_time  [nParticles] = EC_Time[i];
                    ec_path  [nParticles] = EC_Path[i];
                    ec_in    [nParticles] = EC_in  [i];
                    ec_out   [nParticles] = EC_out [i];
                    ec_tot   [nParticles] = EC_tot [i];
                    Mass     [nParticles] = mass   [i];
                    ec_x     [nParticles] = EC_X   [i];
                    ec_y     [nParticles] = EC_Y   [i];
                    ec_z     [nParticles] = EC_Z   [i];
                    ec_u     [nParticles] = EC_U   [i];
                    ec_v     [nParticles] = EC_V   [i];
                    ec_w     [nParticles] = EC_W   [i];
                    charge   [nParticles] = Charge [i];
                    beta     [nParticles] = Beta   [i];
                    theta[nParticles] = Theta[i];
                    phi[nParticles] = Phi[i];
		   
                    // Fill diagnostic histograms for passing protons
                    h2_p_deltaTmom2-> Fill(delta_t   ,mom                   [i]);
                    h2_p_phiTheta2 -> Fill(Phi    [i],Theta                 [i]);
                    h2_p_vzVzCor   -> Fill(targetZ[i],p_vz_corrected-targetZ[i]);
                    h2_p_p_momCor0 -> Fill(mom    [i],mom [i]-p_mom_corrected  );
                    h2_p_p_momCor1 -> Fill(mom    [i],p_mom_corrected/mom   [i]);
                    h2_p_th_pCor0  -> Fill(Theta  [i],mom [i]-p_mom_corrected  );
                    h2_p_th_pCor1  -> Fill(Theta  [i],p_mom_corrected/mom   [i]);
                    h2_p_th_p_cor0 -> Fill(Theta  [i],mom[i],mom[i]-p_mom_corrected);
                    h2_p_th_p_cor1 -> Fill(Theta  [i],mom[i],p_mom_corrected/mom[i]);
                    h2_p_phiVz0    -> Fill(Phi    [i],targetZ               [i]);
                    h2_p_phiVz     -> Fill(Phi    [i],p_vz_corrected           );
                    h2_p_thetaVz0  -> Fill(Theta  [i],targetZ               [i]);
                    h2_p_thetaVz   -> Fill(Theta  [i],p_vz_corrected	   );
                    h2_p_pBeta     -> Fill(mom    [i],Beta                  [i]);
		    
                    nProtons++;
                    nParticles++;
                  }
              }
              // --------------------------------------------------------------------
              // Look specifically for pions +
              else if( (max_pip >= min_pip) && // if we care about pi+
		       (id_guess[i] == 211)&&  // and we have a candidate pi+
		       (fid_params.in_pip_deltaT(pip_delta_t, mom[i], pipdeltat_sig_cutrange)) // Pi+ PID
		       ){
                h2_pip_deltaTmom1 -> Fill(pip_delta_t,mom  [i]);

                // Passed PID for pi+, so check fiducials
                if (!do_fid_pip || (do_fid_pip && fid_params.pFiducialCut(T3_p_mom)))
                  {

                    Part_type[nParticles] = 211;
                    e_deltat [nParticles] = pip_delta_t;
                    mom_x    [nParticles] = T3_p_mom.X();
                    mom_y    [nParticles] = T3_p_mom.Y();
                    mom_z    [nParticles] = T3_p_mom.Z();
                    vtx_z_unc[nParticles] = targetZ  [i];
                    vtx_z_cor[nParticles] = p_vz_corrected;
                    stat_sc  [nParticles] = StatSC [i];
                    stat_dc  [nParticles] = StatDC [i];
                    stat_ec  [nParticles] = StatEC [i];
                    sc_time  [nParticles] = SC_Time[i];
                    sc_path  [nParticles] = SC_Path[i];
                    sc_pad  [nParticles] = SC_Pad[i];
                    ec_time  [nParticles] = EC_Time[i];
                    ec_path  [nParticles] = EC_Path[i];
                    ec_in    [nParticles] = EC_in  [i];
                    ec_out   [nParticles] = EC_out [i];
                    ec_tot   [nParticles] = EC_tot [i];
                    Mass     [nParticles] = mass   [i];
                    ec_x     [nParticles] = EC_X   [i];
                    ec_y     [nParticles] = EC_Y   [i];
                    ec_z     [nParticles] = EC_Z   [i];
                    ec_u     [nParticles] = EC_U   [i];
                    ec_v     [nParticles] = EC_V   [i];
                    ec_w     [nParticles] = EC_W   [i];
                    charge   [nParticles] = Charge [i];
                    beta     [nParticles] = Beta   [i];
                    theta[nParticles] = Theta[i];
                    phi[nParticles] = Phi[i];
		      
                    // Diagnostic histograms
                    h2_pip_pBeta      -> Fill(mom[i]          ,Beta [i]);	
                    h2_pip_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);
		      
                    nPiplus++;
                    nParticles++;
                  }
              }
	      
              // --------------------------------------------------------------------
            }
          // ------------------------------------------------------------------------------------------
          // Test if negative particle
          else if(        (StatSC[i] > 0) &&              // SC status is good for the positive candidate
                          (StatDC[i] > 0) &&              // DC status is good for the positive candidate
                          (Stat  [i] > 0) &&              // Global status is good for the positive candidate
                          (Charge[i] < 0)                 // Charge is negative
                          )
            {
              h2_pim_deltaTmom0 -> Fill(pip_delta_t,mom  [i]);

	      //cout << "outside" << endl;
              if((max_pim >= min_pim) && // if we care about pi-
                 (id_guess[i] == -211)&& // and we have a pi- candidate
                 (fid_params.in_pim_deltaT(pip_delta_t, mom[i], pipdeltat_sig_cutrange)) // Pi- PID
                 ){
		h2_pim_deltaTmom1 -> Fill(pip_delta_t,mom  [i]);
		//cout << "in PID" << endl;
                // No fiducial cuts in place just yet
                if (!do_fid_pim || (do_fid_pim && fid_params.pim_inFidRegion(T3_p_mom)))
		  {		  
		    //cout << "in FID" << endl;
		    Part_type[nParticles] = -211;
		    e_deltat [nParticles] = pip_delta_t;
		    mom_x    [nParticles] = T3_p_mom.X();
		    mom_y    [nParticles] = T3_p_mom.Y();
		    mom_z    [nParticles] = T3_p_mom.Z();
		    vtx_z_unc[nParticles] = targetZ  [i];
		    vtx_z_cor[nParticles] = p_vz_corrected;
		    stat_sc  [nParticles] = StatSC [i];
		    stat_dc  [nParticles] = StatDC [i];
		    stat_ec  [nParticles] = StatEC [i];
		    sc_time  [nParticles] = SC_Time[i];
		    sc_path  [nParticles] = SC_Path[i];
		    sc_pad  [nParticles] = SC_Pad[i];
		    ec_time  [nParticles] = EC_Time[i];
		    ec_path  [nParticles] = EC_Path[i];
		    ec_in    [nParticles] = EC_in  [i];
		    ec_out   [nParticles] = EC_out [i];
		    ec_tot   [nParticles] = EC_tot [i];
		    Mass     [nParticles] = mass   [i];
		    ec_x     [nParticles] = EC_X   [i];
		    ec_y     [nParticles] = EC_Y   [i];
		    ec_z     [nParticles] = EC_Z   [i];
		    ec_u     [nParticles] = EC_U   [i];
		    ec_v     [nParticles] = EC_V   [i];
		    ec_w     [nParticles] = EC_W   [i];
		    charge   [nParticles] = Charge [i];
		    beta     [nParticles] = Beta   [i];
		    theta[nParticles] = Theta[i];
		    phi[nParticles] = Phi[i];
		 
		    h2_pim_pBeta      -> Fill(mom[i]          ,Beta [i]);
		    h2_pim_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);

		    nPiminus++;
		    nParticles++;
		  }
	      }
            }

          // ------------------------------------------------------------------------------------------
          // Test if neutral particle
          else if(        StatEC[i] > 0 && // EC status is good
                          StatDC[i] <=0 && // DC status is not negative
                          StatSC[i] <=0 && // SC status is not negative
                          Stat  [i] > 0 && // Global status is good
                          Charge[i] ==0 && // Charge is neutral
                          Beta  [i] < 1
                          )
            {
	
              EC_Path_corr = fid_params.corrected_path_length( EC_Path[i] , EC_in[i] , EC_out[i] );
              Beta_corr = EC_Path_corr / (c_cm_ns*(EC_Time[i]-t0));
	      
              h2_n_phiTheta0 -> Fill(Phi [i],Theta[i] );
              h2_neu_pBeta   -> Fill(mom [i],Beta_corr);
              h2_n_phiTheta1 -> Fill(Phi [i],Theta[i] );
	
              // --------------------------------------------------------------------
              // Look specifically for neutrons 
              if( (max_n >= min_n) && // if we care about neutrons
		  (Beta_corr < 0.95 ))
                // Don't use: id_guess[i] == 2112 since there is no track.
                // Instead see if ECal hit is neutron like.
                {
                  h1_u_0 -> Fill(EC_U[i]);
                  h1_v_0 -> Fill(EC_V[i]);
                  h1_w_0 -> Fill(EC_W[i]);

                  h1_x_0 -> Fill(EC_X[i]);
                  h1_y_0 -> Fill(EC_Y[i]);
                  h1_z_0 -> Fill(EC_Z[i]);

                  h2_n_ECxy_0 -> Fill(EC_X[i],EC_Y[i]);

                  TVector3 n_ec_xyz(EC_X[i],EC_Y[i],EC_Z[i]);

                  if ((!do_fid_n) || ((do_fid_n) && (fid_params.CutUVW( n_ec_xyz ,10.))))
                    {
                      n_p  = Beta_corr*mN/sqrt(1-Beta_corr*Beta_corr);
		      
                      n_px = n_p*u1.X();
                      n_py = n_p*u1.Y();
                      n_pz = n_p*u1.Z();
		      
                      Part_type[nParticles] = 2112;
                      e_deltat [nParticles] = pip_delta_t;
                      mom_x    [nParticles] = n_px;
                      mom_y    [nParticles] = n_py;
                      mom_z    [nParticles] = n_pz;
                      //mom_x    [nParticles] = T3_p_mom.X();
                      //mom_y    [nParticles] = T3_p_mom.Y();
                      //mom_z    [nParticles] = T3_p_mom.Z();
                      vtx_z_unc[nParticles] = targetZ  [i];
                      vtx_z_cor[nParticles] = p_vz_corrected;
                      stat_sc  [nParticles] = StatSC [i];
                      stat_dc  [nParticles] = StatDC [i];
                      stat_ec  [nParticles] = StatEC [i];
                      sc_time  [nParticles] = SC_Time[i];
                      sc_path  [nParticles] = SC_Path[i];
                      sc_pad  [nParticles] = SC_Pad[i];
                      ec_time  [nParticles] = EC_Time[i];
                      ec_path  [nParticles] = EC_Path[i];
                      ec_in    [nParticles] = EC_in  [i];
                      ec_out   [nParticles] = EC_out [i];
                      ec_tot   [nParticles] = EC_tot [i];
                      Mass     [nParticles] = mass   [i];
                      ec_x     [nParticles] = EC_X   [i];
                      ec_y     [nParticles] = EC_Y   [i];
                      ec_z     [nParticles] = EC_Z   [i];
                      ec_u     [nParticles] = EC_U   [i];
                      ec_v     [nParticles] = EC_V   [i];
                      ec_w     [nParticles] = EC_W   [i];
                      charge   [nParticles] = Charge [i];
                      beta     [nParticles] = Beta_corr ;
                      theta[nParticles] = Theta[i];
                      phi[nParticles] = Phi[i];

                      // Diagnostic histograms
                      h1_u_1 -> Fill(EC_U[i]);
                      h1_v_1 -> Fill(EC_V[i]);
                      h1_w_1 -> Fill(EC_W[i]);
                      h1_x_1 -> Fill(EC_X[i]);
                      h1_y_1 -> Fill(EC_Y[i]);
                      h1_z_1 -> Fill(EC_Z[i]);
                      h2_n_ECxy_1 -> Fill(EC_X[i],EC_Y[i]);
                      h2_n_phiTheta2 -> Fill(Phi[i],Theta[i] );
                      h2_n_pBeta     -> Fill(n_p   ,Beta_corr);
                      nNeutrons++;
                      nParticles++;
                    }
		  
                }

            }

        }

      // --------------------------------------------------------------------------------------------------
      // Fill the output tree
      bool protonsOk = ((max_p < min_p) || ((nProtons >= min_p) && (nProtons <= max_p)));
      bool neutronsOk = ((max_n < min_n) || ((nProtons >= min_n) && (nProtons <= max_n)));
      bool piPlusOk = ((max_pip < min_pip) || ((nProtons >= min_pip) && (nProtons <= max_pip)));
      bool piMinusOk = ((max_pim < min_pim) || ((nProtons >= min_pim) && (nProtons <= max_pim)));

      if( protonsOk && neutronsOk && piPlusOk && piMinusOk)
        outtree->Fill();
    }
  
  cerr << "Finished with the event loop...\n";

  // --------------------------------------------------------------------------------------------------
  // Write the output file
  outfile->cd();
  outtree->Write();

  // Electron diagnostics
  dir_e->cd();
  h1_e_EC_in0       ->Write();
  h1_e_EC_in1       ->Write();
  h1_e_EC_in2       ->Write();
  h1_e_EC_out0      ->Write();
  h1_e_EC_out1      ->Write();
  h1_e_EC_out2      ->Write();
  h1_e_EC_tot0      ->Write();
  h1_e_EC_tot1      ->Write();
  h1_e_EC_tot2      ->Write();
  h2_e_thetaMom0    ->Write();
  h2_e_thetaMom1    ->Write();
  h2_e_thetaMom2    ->Write();
  h1_e_momCor       ->Write();
  h1_e_momCor1      ->Write();
  h2_e_momMomCor    ->Write();
  h2_e_momMomCor1   ->Write();
  h2_e_vzVzCor      ->Write();
  h2_e_Ein_Eout0    ->Write();
  h2_e_Ein_Eout1    ->Write();
  h2_e_Ein_Eout2    ->Write();
  h2_e_xyEC_hit0    ->Write();
  h2_e_xyEC_hit1    ->Write();
  h2_e_xyEC_hit2    ->Write();
  h2_e_p_Etot0      ->Write();
  h2_e_p_Etot1      ->Write();
  h2_e_p_Etot2      ->Write();
  h2_e_p_E0         ->Write();
  h2_e_p_E1         ->Write();
  h2_e_p_E2         ->Write();
  h1_e_Nphe0        ->Write();
  h1_e_Nphe1        ->Write();
  h1_e_Nphe2        ->Write();
  h2_e_phiTheta0    ->Write();
  h2_e_phiTheta1    ->Write();
  h2_e_phiTheta2    ->Write();
  h1_e_vz0          ->Write();
  h1_e_vz           ->Write();
  h1_e_vz_sec10     ->Write();
  h1_e_vz_sec1      ->Write();
  h1_e_vz_sec20     ->Write();
  h1_e_vz_sec2      ->Write();
  h1_e_vz_sec30     ->Write();
  h1_e_vz_sec3      ->Write();
  h1_e_vz_sec40     ->Write();
  h1_e_vz_sec4      ->Write();
  h1_e_vz_sec53     ->Write();
  h1_e_vz_sec5      ->Write();
  h1_e_vz_sec60     ->Write();
  h1_e_vz_sec6      ->Write();
  h2_e_phiVz0       ->Write();
  h2_e_phiVz        ->Write();
  h2_e_thetaVz0     ->Write();
  h2_e_thetaVz      ->Write();
  // ---
  dir_pos->cd();
  h2_pos_pBeta      ->Write();
  h2_p_pBeta        ->Write();
  h2_pip_pBeta      ->Write();
  // ---
  dir_p->cd();
  h1_p_mass         ->Write();
  h2_p_pMass        ->Write();	
  h2_p_p_momCor0    ->Write();
  h2_p_vzVzCor      ->Write();
  h2_p_phiVz0       ->Write(); 
  h2_p_phiVz        ->Write();
  h2_p_thetaVz0     ->Write();
  h2_p_thetaVz      ->Write();
  h2_p_deltaTmom0   ->Write();
  h2_p_deltaTmom1   ->Write();
  h2_p_deltaTmom2   ->Write();
  h2_p_phiTheta0    ->Write();
  h2_p_phiTheta1    ->Write();
  h2_p_phiTheta2    ->Write();
  // ---
  dir_pip->cd();
  h2_pip_deltaTmom0 ->Write();
  h2_pip_deltaTmom1 ->Write();
  h2_pip_deltaTmom2 ->Write();
  // ---
  dir_neut->cd();
  h2_neu_pBeta -> Write();
  // ---
  dir_n->cd();
  h2_n_phiTheta0->Write();
  h2_n_phiTheta1->Write();
  h2_n_phiTheta2->Write();
  h2_n_pBeta->Write();
  h1_u_0->Write();
  h1_v_0->Write();
  h1_w_0->Write();
  h1_u_1->Write();
  h1_v_1->Write();
  h1_w_1->Write();
  h1_x_0->Write();
  h1_y_0->Write();
  h1_z_0->Write();
  h1_x_1->Write();
  h1_y_1->Write();
  h1_z_1->Write();
  h2_n_ECxy_0->Write();
  h2_n_ECxy_1->Write();
  // ---
  dir_pim->cd();  
  h2_pim_pBeta->Write();
  h2_pim_deltaTmom0->Write();
  h2_pim_deltaTmom1->Write();
  h2_pim_deltaTmom2->Write();

  // Clean up
  outfile->Close();

  return 0;
}






// Rey's artistic code. We'll package this into another program at some point
/*
// --------------------------------------------------------------------------------------------------
// Editing histograms
h1_e_Nphe1  -> SetLineColor(2);
h1_e_Nphe2  -> SetLineColor(8);
h1_e_EC_in1 -> SetLineColor(2);
h1_e_EC_in2 -> SetLineColor(8);
h1_e_EC_out1-> SetLineColor(2);
h1_e_EC_out2-> SetLineColor(8);
h1_e_EC_tot1-> SetLineColor(2);
h1_e_EC_tot2-> SetLineColor(8);
h1_Xb1      -> SetLineColor(2);
h1_e_vz     -> SetLineColor(2);
// ---
TLegend * leg = new TLegend(0.5,0.5,0.8,0.8);
leg -> AddEntry(h1_e_Nphe0,"before cuts"       );
leg -> AddEntry(h1_e_Nphe1,"after PID cuts"    );
leg -> AddEntry(h1_e_Nphe2,"after PID+Fiducial");
// --------------------------------------------------------------------------------------------------
TCanvas *c1 = new TCanvas("c1");
h1_e_Nphe0 -> Draw();
h1_e_Nphe1 -> Draw("same");
h1_e_Nphe2 -> Draw("same");
leg          -> Draw("same");

TCanvas *c2 = new TCanvas("c2");
h1_e_EC_in0 -> Draw();
h1_e_EC_in1 -> Draw("same");
h1_e_EC_in2 -> Draw("same");
leg           -> Draw("same");

TCanvas *c3 = new TCanvas("c3");
h1_e_EC_out0 -> Draw();
h1_e_EC_out1 -> Draw("same");
h1_e_EC_out2 -> Draw("same");
leg            -> Draw("same");

TCanvas *c4 = new TCanvas("c4");
h1_e_EC_tot0 -> Draw();
h1_e_EC_tot1 -> Draw("same");
h1_e_EC_tot2 -> Draw("same");
leg            -> Draw("same");

TCanvas *c5 = new TCanvas("c5");
c5 -> Divide(3,1);
c5 -> cd(1);	h2_e_thetaMom0 -> Draw("COLZ"); 
c5 -> cd(2);	h2_e_thetaMom1 -> Draw("COLZ");
c5 -> cd(3);	h2_e_thetaMom2 -> Draw("COLZ");

TCanvas *c6 = new TCanvas("c6");
h2_e_thetaMom3 -> Draw("COLZ");

TCanvas *c7 = new TCanvas("c7");
c7 -> Divide(2,1);
c7 -> cd(1);	h1_e_momCor -> Draw();
c7 -> cd(2);    h1_e_momCor1-> Draw();

TCanvas *c8 = new TCanvas("c8");
c8 -> Divide(2,1);
c8 -> cd(1);    h2_e_momMomCor -> Draw("COLZ");
c8 -> cd(2);    h2_e_momMomCor1-> Draw("COLZ");

TCanvas *c9 = new TCanvas("c9");
c9 -> Divide(3,2);
c9 -> cd(1);	h1_e_momCor_sec1 -> Draw();
c9 -> cd(2);	h1_e_momCor_sec2 -> Draw();
c9 -> cd(3);	h1_e_momCor_sec3 -> Draw();
c9 -> cd(4);	h1_e_momCor_sec4 -> Draw();
c9 -> cd(5);	h1_e_momCor_sec5 -> Draw();
c9 -> cd(6);	h1_e_momCor_sec6 -> Draw();

TCanvas *c10 = new TCanvas("c10");
c10 -> Divide(3,2);
c10 -> cd(1);	h2_e_momMomCor_sec1 -> Draw("COLZ");
c10 -> cd(2);	h2_e_momMomCor_sec2 -> Draw("COLZ");
c10 -> cd(3);	h2_e_momMomCor_sec3 -> Draw("COLZ");
c10 -> cd(4);	h2_e_momMomCor_sec4 -> Draw("COLZ");
c10 -> cd(5);	h2_e_momMomCor_sec5 -> Draw("COLZ");
c10 -> cd(6);	h2_e_momMomCor_sec6 -> Draw("COLZ");

TCanvas *c11 = new TCanvas("c11");
h2_e_vzVzCor -> Draw("COLZ");

TCanvas *c12 = new TCanvas("c12");
h1_Xb0 -> Draw();
h1_Xb1 -> Draw("same");

TCanvas *c13 = new TCanvas("c13");
c13 -> Divide(3,1);
c13 -> cd(1);	h2_e_Ein_Eout0 -> Draw("COLZ");
c13 -> cd(2);	h2_e_Ein_Eout1 -> Draw("COLZ");
c13 -> cd(3);	h2_e_Ein_Eout2 -> Draw("COLZ");

TCanvas *c14 = new TCanvas("c14");
c14 -> Divide(3,1);
c14 -> cd(1);    h2_e_Ein_Eout_0 -> Draw("COLZ");
c14 -> cd(2);    h2_e_Ein_Eout_1 -> Draw("COLZ");
c14 -> cd(3);    h2_e_Ein_Eout_2 -> Draw("COLZ");

TCanvas *c15 = new TCanvas("c15");
c15 -> Divide(3,1);
c15 -> cd(1);	h2_e_xyEC_hit0 -> Draw("COLZ");
c15 -> cd(2);	h2_e_xyEC_hit1 -> Draw("COLZ");
c15 -> cd(3);	h2_e_xyEC_hit2 -> Draw("COLZ");

TCanvas *c16 = new TCanvas("c16");
c16 -> Divide(3,1);
c16 -> cd(1);	h2_e_p_Etot0 -> Draw("COLZ");
c16 -> cd(2);	h2_e_p_Etot1 -> Draw("COLZ");
c16 -> cd(3);	h2_e_p_Etot2 -> Draw("COLZ");

TCanvas *c17 = new TCanvas("c17");
c17 -> Divide(3,1);
c17 -> cd(1);	h2_e_p_E0 -> Draw("COLZ");
c17 -> cd(2);	h2_e_p_E1 -> Draw("COLZ");
c17 -> cd(3);	h2_e_p_E2 -> Draw("COLZ");

TCanvas *c18 = new TCanvas("c18");	
h2_e_phiTheta0 -> Draw("COLZ");
TCanvas *c19 = new TCanvas("c19");
h2_e_phiTheta1 -> Draw("COLZ");
TCanvas *c20 = new TCanvas("c20");
h2_e_phiTheta2 -> Draw("COLZ");

TCanvas *c21 = new TCanvas("c21");
h1_e_vz0 -> Draw();
h1_e_vz  -> Draw("same");

TCanvas *c22 = new TCanvas("c22");
c22 -> Divide(3,2);
c22 -> cd(1);	h1_e_vz_sec10 -> Draw();
c22 -> cd(2);	h1_e_vz_sec20 -> Draw();
c22 -> cd(3);	h1_e_vz_sec30 -> Draw();
c22 -> cd(4);	h1_e_vz_sec40 -> Draw();
c22 -> cd(5);	h1_e_vz_sec53 -> Draw();
c22 -> cd(6);	h1_e_vz_sec60 -> Draw();

TCanvas *c23 = new TCanvas("c23");
c23 -> Divide(3,2);
c23 -> cd(1);	h1_e_vz_sec1 -> Draw();
c23 -> cd(2);	h1_e_vz_sec2 -> Draw();
c23 -> cd(3);	h1_e_vz_sec3 -> Draw();
c23 -> cd(4);	h1_e_vz_sec4 -> Draw();
c23 -> cd(5);	h1_e_vz_sec5 -> Draw();
c23 -> cd(6);	h1_e_vz_sec6 -> Draw();

TCanvas *c24 = new TCanvas("c24");
c24 -> Divide(2,1);
c24 -> cd(1);	h2_e_phiVz0 -> Draw("COLZ");
c24 -> cd(2);	h2_e_phiVz  -> Draw("COLZ");

TCanvas *c25 = new TCanvas("c25");
c25 -> Divide(2,1);
c25 -> cd(1);	h2_e_thetaVz0 -> Draw("COLZ");
c25 -> cd(2);	h2_e_thetaVz  -> Draw("COLZ");

TCanvas *c26 = new TCanvas("c26");
h2_pos_pBeta -> Draw("COLZ");

TCanvas *c27 = new TCanvas("c27");
h2_p_pBeta -> Draw("COLZ");

TCanvas *c28 = new TCanvas("c28");
c28 -> Divide(2,1);
c28 -> cd(1);	h1_p_mass ->Draw();
c28 -> cd(2);	h2_p_pMass->Draw("COLZ");

TCanvas *c29 = new TCanvas("c29");
c29 -> Divide(2,1);
c29 -> cd(1);	h2_p_p_momCor0-> Draw("COLZ");
c29 -> cd(2);   h2_p_p_momCor1-> Draw("COLZ");

TCanvas *c30 = new TCanvas("c30");
c30 -> Divide(2,1);
c30 -> cd(1);	h2_p_th_pCor0-> Draw("COLZ");
c30 -> cd(2);	h2_p_th_pCor1-> Draw("COLZ");

TCanvas *c31 = new TCanvas("c31");
c31 -> Divide(2,1);
c31 -> cd(1);	h2_p_th_p_cor0 -> Draw("COLZ");
c31 -> cd(2);	h2_p_th_p_cor1 -> Draw("COLZ");

TCanvas *c32 = new TCanvas("c32");
h2_p_vzVzCor -> Draw("COLZ");

TCanvas *c33 = new TCanvas("c33");
c33 -> Divide(2,1);
c33 -> cd(1);	h2_p_phiVz0 -> Draw("COLZ");
c33 -> cd(2);	h2_p_phiVz  -> Draw("COLZ");

TCanvas *c34 = new TCanvas("c34");
c34 -> Divide(2,1);
c34 -> cd(1);	h2_p_thetaVz0 -> Draw("COLZ");
c34 -> cd(2);	h2_p_thetaVz  -> Draw("COLZ");

TCanvas *c35 = new TCanvas("c35");
h2_p_deltaTmom0 -> Draw("COLZ");

TCanvas *c36 = new TCanvas("c36");
h2_p_deltaTmom1 -> Draw("COLZ");

TCanvas *c37 = new TCanvas("c37");
h2_p_deltaTmom2 -> Draw("COLZ");

TCanvas *c38 = new TCanvas("c38");
h2_p_phiTheta0 -> Draw("COLZ");

TCanvas *c39 = new TCanvas("c39");
h2_p_phiTheta1 -> Draw("COLZ");

TCanvas *c40 = new TCanvas("c40");
h2_p_phiTheta2 -> Draw("COLZ");

TCanvas *c41 = new TCanvas("c41");
h2_pip_deltaTmom0 -> Draw("COLZ");

TCanvas *c42 = new TCanvas("c42");
h2_pip_deltaTmom1 -> Draw("COLZ");

TCanvas *c43 = new TCanvas("c43");
h2_pip_deltaTmom2 -> Draw("COLZ");

TCanvas *c44 = new TCanvas("c44");
h2_neu_pBeta -> Draw("COLZ");

TCanvas *c45 = new TCanvas("c45");
h2_n_pBeta -> Draw("COLZ");

TCanvas *c46 = new TCanvas("c46");
h2_n_phiTheta0 -> Draw("COLZ");

TCanvas *c47 = new TCanvas("c47");
h2_n_phiTheta1 -> Draw("COLZ");

TCanvas *c48 = new TCanvas("c48");
h2_n_phiTheta2 -> Draw("COLZ");

TCanvas *c49 = new TCanvas("c49");
c49 -> Divide(3,2);
c49 -> cd(1);	h1_u_0 -> Draw();
c49 -> cd(2);   h1_v_0 -> Draw();
c49 -> cd(3);   h1_w_0 -> Draw();
c49 -> cd(4);   h1_u_1 -> Draw();
c49 -> cd(5);   h1_v_1 -> Draw();
c49 -> cd(6);   h1_w_1 -> Draw();

TCanvas *c50 = new TCanvas("c50");
c50 -> Divide(3,2);
c50 -> cd(1);   h1_x_0 -> Draw();
c50 -> cd(2);   h1_y_0 -> Draw();
c50 -> cd(3);   h1_z_0 -> Draw();
c50 -> cd(4);   h1_x_1 -> Draw();
c50 -> cd(5);   h1_y_1 -> Draw();
c50 -> cd(6);   h1_z_1 -> Draw();

TCanvas *c51 = new TCanvas("c51");
h2_n_ECxy_0 -> Draw("COLZ");

TCanvas *c52 = new TCanvas("c52");
h2_n_ECxy_1 -> Draw("COLZ");

TCanvas *c53 = new TCanvas("c53");
h2_pim_deltaTmom0 -> Draw("COLZ");

TCanvas *c54 = new TCanvas("c54");
h2_pim_deltaTmom1 -> Draw("COLZ");

TCanvas *c55 = new TCanvas("c55");
h2_pim_deltaTmom2 -> Draw("COLZ");

// --------------------------------------------------------------------------------------------------
// Print histograms on a pdf file

c1  -> Print(Form("../plots/plots_%d.pdf(",tab_run),"pdf");
c2  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c3  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c4  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c5  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c6  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c7  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c8  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c9  -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c10 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c11 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c12 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c13 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c14 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c15 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c16 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c17 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c18 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c19 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c20 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c21 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c22 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c23 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c24 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c25 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c26 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c27 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c28 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c29 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c30 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c31 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c32 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c33 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c34 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c35 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c36 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c37 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c38 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c39 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c40 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c41 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");	
c42 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c43 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c44 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c45 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c46 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c47 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c48 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c49 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c50 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c51 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c52 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c53 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c54 -> Print(Form("../plots/plots_%d.pdf",tab_run) ,"pdf");
c55 -> Print(Form("../plots/plots_%d.pdf)",tab_run),"pdf");

*/
