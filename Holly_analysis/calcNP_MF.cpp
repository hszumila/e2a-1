#include "TH1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TROOT.h"
#include <vector>
#include "unistd.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TCut.h"
#include "Riostream.h"
#include "TNtupleD.h"
#include "TCut.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Fiducial.h"
#include "Acceptance.h"


//define some constants
const double amu     = 0.931494;
const double Me      = .000511;         	// Electron mass in GeV
const double Mp      = 0.93827;         	// Proton   mass in GeV
const double Mn      = 0.93957;         	// Neutron  mass in GeV
const double BE_2H   = 0.00222;			// 2H  binding energy in GeV
const double BE_3H   = 0.00848;			// 3H  binding energy in GeV
const double BE_3He  = 0.00772;                 // 3He binding energy in GeV
const double BE_4He  = 0.02830;                 // 4He binding energy in GeV
const double BE_11B  = 0.07620;			// 11B binding energy in GeV
const double BE_11C  = 0.07344;			// 11C binding energy in GeV
const double BE_12C  = 0.09215;                 // 12C binding energy in GeV
const double M2H     =   Mp +   Mn -BE_2H ;
const double M3H     =   Mp + 2*Mn -BE_3H ;
const double M3He    = 2* Mp+Mn -BE_3He;        // 3He mass in GeV
const double M4He    = 2*(Mp+Mn)-BE_4He;        // 4He mass in GeV
const double M11B    = 5*Mp + 6*Mn -BE_11B;
const double M11C    = 6*Mp + 5*Mn -BE_11C;
const double M12C    = 6*(Mp+Mn)-BE_12C;        // 12C mass in GeV
const double pi      = 3.14159;
const double cc      = 2.99792458E+8;

//vertex cuts
const double vtx_min_C = 4.0;//cm
const double vtx_max_C = 7.0;//cm
const double vtx_min_He = -2.5;//cm
const double vtx_max_He = -0.5;//cm

//SRC cuts defined:
const double xbP_cut = 1.2;
const double xbN_cut = 1.1;
const double thetapqSRC_cut = 25.0;//deg
const double pq_min = 0.62;
const double pq_max_n = 1.1;
const double pq_max_p = 0.96;
const double mmiss_nomCut_n = 1.15;//1.1;//1.15;//1.1;//1.2;//1.15;//1.175;
const double pmissSRC_nomCut_n = 0.36;//0.3;//0.36;//0.3;//0.38;//0.36;//0.38;
const double mmiss_nomCut_p = 1.1;
const double pmissSRC_nomCut_p = 0.3;

//MF cuts defined:
const double y_cut_min = -0.05;
const double y_cut_max = 0.2;
const double nu_min = 0.9;
const double nu_max = 1.6;
const double thetapqMF_cut = 7.0;//deg
const double emiss_nomCut = 0.08;
const double pmissMF_nomCut = 0.22;
const double emiss_nomCut_n = 0.22;
const double pmissMF_nomCut_n = 0.28;


const double Q2_min = 0.0;//1.5;



//// Q2 in GeV^2
////
// The parameterization formula returns the uncertainty devided by G(0)*GD, where 
//  GD(Q2) = 1./(1+Q2/0.71)^2
// and GEp(0) = 1, GMp(0) = 2.79284356, GEn(0) = 1, GMn(0) = -1.91304272,
//
// The parameterization formula for the Form Factor value is:
//  $$ GN(z) = sum_{i=0}^{N=12}(a_i * z^i) 
// Note that the return value has been divided by (G(Q2=0)*G_Dip)
//
// The parameterization formula for the Form Factor error is:
// $$ log_{10}\frac{\delta G}{G_D} = (L+c_0)\Theta_a(L_1-L) 
//                                 +\sum_{i=1}^{N}(c_i+d_i L)[\Theta_a(L_i-L)-\Theta_a(L_{i+1}-L)]
//                                 +log_{10}(E_{\inf})\Theta_a(L-L_{N+1})$$
// where $L=log_{10}(Q^2)$, $\Theta_{a}(x)=[1+10^{-ax}]^{-1}$. $a=1$.

double GetFF(const int kID, const double kQ2, int rval=1){//, double *GNGD_Fit, double* GNGD_Err){/*{{{*/
    // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  //rval: 1-FF, 2-error

  double GNGD_Fit, GNGD_Err;
    if (kID<1 || kID>4){
        cerr<<"*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]"<<endl;
        //GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
        GNGD_Fit = -1000;  GNGD_Err = -1000;
        return -1;
    }
    ////////////////////////////////////////////////
    //// z-Expansion Parameters for Form Factor Values
    /////////////////////////////////////////////////*{{{*/
    const double GN_Coef_Fit[4][13] ={
        {0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, -2.28689474187,  1.12663298498, 1.25061984354,-3.63102047159, 4.08221702379,  0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103}, /*GEp*/
        {0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, -1.35641843888, 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    8.70740306757, -5.70699994375, 1.28081437589}, /*GMp*/
        {0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0}, /*GEn*/
        {0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0} /*GMn*/
    };/*}}}*/

    ////////////////////////////////////////////////
    //// Parameters for Form Factor Errors
    ////////////////////////////////////////////////  
    const double parL[4][2] ={
        {-0.97775297,  0.99685273}, //GEp
        {-0.68452707,  0.99709151}, //GMp
        {-2.02311829, 1.00066282}, //GEn
        {-0.20765505, 0.99767103}, //GMn
    };
    const double parM[4][15] = {
        {  -1.97750308e+00,  -4.46566998e-01,   2.94508717e-01,   1.54467525e+00,
            9.05268347e-01,  -6.00008111e-01,  -1.10732394e+00,  -9.85982716e-02,
            4.63035988e-01,   1.37729116e-01,  -7.82991627e-02,  -3.63056932e-02,
            2.64219326e-03,   3.13261383e-03,   3.89593858e-04}, //GEp

        {  -1.76549673e+00,   1.67218457e-01,  -1.20542733e+00,  -4.72244127e-01,
            1.41548871e+00,   6.61320779e-01,  -8.16422909e-01,  -3.73804477e-01,
            2.62223992e-01,   1.28886639e-01,  -3.90901510e-02,  -2.44995181e-02,
            8.34270064e-04,   1.88226433e-03,   2.43073327e-04}, //GMp

        {  -2.07343771e+00,   1.13218347e+00,   1.03946682e+00,  -2.79708561e-01,
           -3.39166129e-01,   1.98498974e-01,  -1.45403679e-01,  -1.21705930e-01,
            1.14234312e-01,   5.69989513e-02,  -2.33664051e-02,  -1.35740738e-02,
            7.84044667e-04,   1.19890550e-03,   1.55012141e-04}, //GEn

        {  -2.07087611e+00,   4.32385770e-02,  -3.28705077e-01,   5.08142662e-01,
            1.89103676e+00,   1.36784324e-01,  -1.47078994e+00,  -3.54336795e-01,
            4.98368396e-01,   1.77178596e-01,  -7.34859451e-02,  -3.72184066e-02,
            1.97024963e-03,   2.88676628e-03,   3.57964735e-04} //GMn
    };
    const double parH[4][3] = {
        {0.78584754,  1.89052183, -0.4104746}, //GEp
        {0.80374002,  1.98005828, -0.69700928}, //GMp
        {0.4553596,  1.95063341,  0.32421279}, //GEn
        {0.50859057, 1.96863291,  0.2321395} //GMn 
    };
    

    //// Apply the z-expansion formula for form factor
    const double tcut = 0.0779191396 ;
    const double t0 = -0.7 ;
    double z = (sqrt(tcut+kQ2)-sqrt(tcut-t0))/(sqrt(tcut+kQ2)+sqrt(tcut-t0)) ;
    double GNQ2 = 0.0;
    for (int i=0;i<13;i++) GNQ2 += GN_Coef_Fit[kID-1][i] * pow(z, i);
    double GDip= pow(1./(1. + kQ2/0.71), 2);
    //GNGD_Fit[0] = GNQ2 / GDip; //Note that Coef_Fit have been divided by mu_p or mu_n
    GNGD_Fit = GNQ2 ;/// GDip; //Note that Coef_Fit have been divided by mu_p or mu_n
    

    //// Apply the parameterization formula for error
    double lnQ2 = log10(kQ2);
    double lnGNGD_Err=0.0;
    if (kQ2<1e-3)
        lnGNGD_Err = parL[kID-1][0] + parL[kID-1][1]*lnQ2;
    else if (kQ2>1e2)
        lnGNGD_Err = parH[kID-1][0]*sqrt(lnQ2-parH[kID-1][1]) + parH[kID-1][2];
    else{
        for (int i=0; i<15;i++) lnGNGD_Err += parM[kID-1][i] * pow(lnQ2,i);
    }
    //GNGD_Err[0] = pow(10.,(lnGNGD_Err));    //LOG10(dG/G(0)/GD);
    GNGD_Err = pow(10.,(lnGNGD_Err));    //LOG10(dG/G(0)/GD);
    
    if (rval==2){
      return GNGD_Err;
    }
    else{return GNGD_Fit;}
}

double y_scale_func(TVector3 y_q, double enu, double y_Mtar, double y_Mres, double y_Mnuc){
	// See equation 9a from "Scaling in inclusive electron-nucleus scattering
	// D.B. Day, J.S. McCarthy, T.W. Donnelly, and I. Sick
	// Annu. Rev. Nucl. Part. Sci 1990. 40: 357-409
	double q3       = y_q.Mag();
	double nu       = enu;
	double W        = sqrt( pow(y_Mtar+nu,2) - pow(q3,2) );
	double Lambda   = ( pow(y_Mres,2) - pow(y_Mnuc,2) + pow(W,2) )/2.;
	return ( (y_Mtar+nu)*sqrt( pow(Lambda,2) - pow(y_Mres*W,2) ) - q3*Lambda )/pow(W,2);     
}

double fn_Mmiss(double omega, double Mnuc, double Enuc, TVector3 q_vector, TVector3 p_vector){
  return sqrt(pow(omega+2*Mnuc-Enuc,2.0) - (q_vector - p_vector).Mag2());
}

double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;								// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;								// Missing energy
}

double f_delta_p(double p){
	// Returns value of delta_p as a function of p.
	// This function is used for neutron momentum resolution.
        // 0.392319c	-0.0967682	0.022168

  //2.2,4.4 GeV 3He:
  double ec_extr_res = 0.091479;//0.098159;//0.091479;//0.098159;//-0.09;//0967682;
  double dt = 0.3292;//0.3532;//0.3292;//0.3532;//0.392319;//0.625                 //ec_time_res;
	const double d      = 516.9152;       //cm, distance from interaction vtx to EC
	const double Mn     = 0.939565378;    // Neutron mass in GeV
	const double c_cm_s = 2.99792458E+10; // cm/s, speed of light
	const double ns_2_s = 1E-09;

	double val = p*p/(Mn*Mn)*sqrt(Mn*Mn+p*p)*dt/d*c_cm_s*ns_2_s;
	return sqrt(val*val/p/p+ec_extr_res*ec_extr_res)*p;
}

float Coulomb_Corr_MeV(int Z, int A){
	// See Effective Momentum Approximation (EMA) formalism
	// https://www.jlab.org/indico/event/206/session/27/contribution/12/material/slides/0.pdf

	const float alpha  = 0.0072973525664;
	const float hbar_c = 197.327; // MeV*fm

	float R0_fm = 1.1*pow((float)A,1/3.) + 0.86*pow((float)A,-1/3.); // nuclear radius [fm]
	float DV = 3*((float)Z-1)*alpha*hbar_c/(2.*R0_fm);
	float DE = 0.775*DV;

	return DE;
}

TVector3 P_Coulomb_Corrected(TVector3 P_Uncorrected, int Z, int A){
	float Correction_GeV = Coulomb_Corr_MeV( Z , A )/1000.;
	double E_uncorr = sqrt(P_Uncorrected.Mag2()+pow(Mp,2.0));
	float E_corr  = E_uncorr - Correction_GeV;
	float New_Mag = sqrt(E_corr*E_corr-Mp*Mp);

	TVector3 unit = P_Uncorrected.Unit();
	unit.SetMag(New_Mag);

	return unit;
}

TVector3 correctPn(TVector3 P, double eIn, double eOut){
  double a,b,c;
  
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    a = 0.06;
    b = -0.17;
    c = 0.08;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    a = 0.02;
    b = -0.05;
    c = 0.05;
  }
  //both
  else {
    a = 0.09;
    b = -0.25;
    c = 0.08;
  }

  double mom = P.Mag();
  double newMom = mom - mom*(a+b*mom+c*mom*mom);
  
  TVector3 unit = P.Unit();
  unit.SetMag(newMom);
  
  return unit;
}


bool p_points_to_EC_fid(TVector3 pm, double vtx_z, double dist){
	double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;

	const double ec_theta =  65.*pi/180.; // degrees
	const double z0       =563.1; // cm

	// Using law of sines to determine distance from vertex to EC
	pm_theta    = pm.Theta();
	third_theta = (pi - ec_theta - pm_theta);

	r = (z0-vtx_z)*sin(ec_theta)/sin(third_theta);

	double Xmiss = r*sin(pm.Theta())*cos(pm.Phi());
	double Ymiss = r*sin(pm.Theta())*sin(pm.Phi());

	pm_phi = 180./pi*pm.Phi();
	if (pm_phi < -30.) pm_phi += 360.;
	int sec = (int)(pm_phi+30)/60;
	if (sec>5) sec = 5;
	if (sec<0) sec = 0;

	angle = pi/180.*60.*(sec);

	x_rot = Xmiss*cos(angle) + Ymiss*sin(angle);
	y_rot = Xmiss*sin(angle) - Ymiss*cos(angle);

	d_intercept = dist/cos(atan(1.73));

	return((x_rot<390.-dist)&&(x_rot>1.73*abs(y_rot)+55.+d_intercept));

}
//NN=1 for proton, 2 for neutron
double getSigmaWeight(double thetaE, double QSq, int NN){
  int ii;
  double partMass;
  if (NN==1){
    ii = 1;
    partMass = Mp;
  }
  else {
    ii = 3;
    partMass = Mn;
  }
  
  double GE = GetFF(ii,QSq,1);
  double GM = GetFF(ii+1,QSq,1);
  double tau = QSq/(4.0*partMass*partMass);
  double epsilon = 1./(1.+2.*(1+tau)*pow(tan(thetaE/2.),2.0));
  double weight = (epsilon/tau)*pow(GE,2.0)+pow(GM,2.0);

  return weight;
}

//get neutron efficiency (from momentum):
double getNEff(double pm){
  double value=0;
  //if (pm<1.7){value = 0.167707*pm -0.0878069;}
  //else {value = 0.167707*1.7 -0.0878069;;}

  //if (pm<1.6){value = 0.32*pm -0.09;}
  //else {value = 0.32*1.6 -0.09;}
  //return value;
  /*
  double par[5] = {1.1,0.162734,-0.0820586,1.74,0.453669};
  if     (pm <  par[0]){value = par[1]*pm    +par[2];}
  else if(pm >= par[0] && pm<par[3]){value = par[4]*pm+(par[1]-par[4])*par[0]+par[2];}
  else {value = par[4]*par[3]+(par[1]-par[4])*par[0]+par[2];}
  return value;
  */
  if (pm<1.65){value = 0.415133*pm -0.188579;}
  else {value = 0.415133*1.65 -0.188579;}
  return value;
}

  


int main(){

  int event_p = 0;
  int event_n = 0;
  int stats_p = 0;
  int stats_n = 0;
  double tgt_min, tgt_max, tgt_M, tgt_R;
  int tgt_Z, tgt_A;
  const int maxParticles=50;
  TFile *fout = new TFile("output/out_np_ratio.root","RECREATE");
  
  TCanvas *canvas = new TCanvas("canvas","output plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLogy();
  std::string pdf_file_name = "output/out_np_ratio.pdf";
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  //input 3He
  TString target = "4He";
  //tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 3; tgt_M = M3He; tgt_R = M2H;
  tgt_min = tgt_min = vtx_min_He; tgt_max = vtx_max_He; tgt_Z = 2; tgt_A = 4; tgt_M = M4He; tgt_R = M3H;
  //tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C; tgt_R = M11B;
  Fiducial fid_params(4461,2250,5996,"4He",true);

  cout<<"target min: "<<tgt_min<<" target max: "<<tgt_max<<endl;

  //open the 3He proton file
  TFile *f_3p = new TFile("../../../data/19jun20/skim_He4_p_Jun20.root");
  TTree * intree = (TTree*)f_3p->Get("T");
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles], ec_in[maxParticles],ec_out[maxParticles];
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

  double theta_prev = 0;

  //get (e,e'p)
  //for (int event =0 ; event < intree->GetEntries() ; event++)
  for (int event =0 ; event < 100000 ; event++)
  //for (int event =0 ; event < 100000 ; event++)
    {
       if (event%10000==0)
      	cout << "Working on event " << event <<endl;

      intree->GetEvent(event);
      if (QSq<Q2_min){continue;}

       // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //mean field stuff:
      double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at protons
	  if (type[part] != 2212)
	    continue;

	  //vertex cut-target dep.
	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Angle between p and q
	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
	  if (theta_pq>thetapqMF_cut){continue;}
	  //if (theta_pq>25.0){continue;}
	  
	  //first proton of interest in event
	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

	  /////////////////////////////////////////////////////////////////////
	  //Here I selected the proton with angle closest to q-vector:	  
	  //current p/q is larger and angle is smaller than previous
	  if (theta_pq<theta_prev){leadingID=part;}
	  /////////////////////////////////////////////////////////////////////	  
	 
	}

      //Now we have a leading proton!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //if (mom_lead.Mag()>2.34){continue;}
	  
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 

	  //apply proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();

	  //check if in good acceptance
	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}


	  //cout<<" y_sc: "<<y_sc<<" Nu: "<<Nu<<"angle: "<<mom_lead.Angle(q)*180./M_PI <<" Pmiss: "<<mom_miss.Mag()<<" Emiss: "<<emiss<<endl;
	  
	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapqMF_cut && mom_miss.Mag()<pmissMF_nomCut && emiss<emiss_nomCut){
	    //cout<<"event survived here"<<endl;

	    double weight = getSigmaWeight(mom_e.Theta(), QSq, 1);
	   
	    Acceptance myMap("He4", 4461, 2250, "p");
	    double p_eff = myMap.get_acc(mom_lead);

	    event_p += 1./p_eff * 1./weight;
	    stats_p += 1;



	  }

	}
    }
  
  intree->Delete();
  
  //get (e,e'n)
  //open the 3He neutron file
  TFile *f_3n = new TFile("../../../data/19jun20/skim_He4_n_Jun20.root");
  intree = (TTree*)f_3n->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  //for (int event =0 ; event < intree->GetEntries() ; event++)
    for (int event =0 ; event < 100000 ; event++)
    {
     
      if (event%10000==0){cout << "Working on event " << event <<endl;}

      intree->GetEvent(event);
      if (QSq<Q2_min){continue;}

       // Get the electron and q vectors
      TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
      double Ebeam = Nu + mom_e.Mag();
      TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
      
      //mean field stuff:
      double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

      // Loop over particles looking for a proton most likely to pass lead cuts
      int leadingID = -1;
      for (int part=1 ; part < nParticles ; part++)
	{   
	  // Only look at protons
	  if (type[part] != 2112)
	    continue;

	  //vertex cut-target dep.
	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
	  
	  // Get a 3-vector for this proton's momentum
	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

	  // Angle between p and q
	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
	  if (theta_pq>thetapqMF_cut){continue;}
	  
	  //first proton of interest in event
	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

	  /////////////////////////////////////////////////////////////////////
	  //Here I selected the proton with angle closest to q-vector:	  
	  //current p/q is larger and angle is smaller than previous
	  if (theta_pq<theta_prev){leadingID=part;}
	  /////////////////////////////////////////////////////////////////////	  
	  	  
	}

      //Now we have a leading neutron!
      if (leadingID > 0)
	{
	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
	  //mom_lead = correctPn(mom_lead,ec_in[leadingID],ec_out[leadingID]);
	  mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

	  //calculate the nominal values
	  TVector3 mom_miss = mom_lead - q;
	  double theta = mom_lead.Theta()*180./M_PI;
	  double phi = mom_lead.Phi();
	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
	  
	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
	  phi *= 180./M_PI; 
	  
	  //apply proton fiducial cuts:
	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}

	  //apply EC 10cm edge projection cut:
	  //this doesn't work, but I don't know why:
	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
	  // Cutting top of EC
	  if(theta>44.) continue;
	  if(mom_lead.Mag()>2.34) continue;
	  
	  //passing lead proton cuts
	  double p_over_q = mom_lead.Mag() / q.Mag();

	  //check if in good acceptance
	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}


	  //select the 3He neutron sample
	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapqMF_cut && mom_miss.Mag()<pmissMF_nomCut_n && emiss<emiss_nomCut_n){
	    double weight = getSigmaWeight(mom_e.Theta(), QSq, 2);	    
	    double n_eff = getNEff(mom_lead.Mag());

	    event_n += 1./(n_eff) * 1./weight;
	    stats_n += 1;


	  }	  
	}//end leading
    }//end event


  double ratio = (double) event_n/event_p;
  double stats = sqrt(pow(1./stats_n,2.0)+pow(1./stats_p,2.0));
  cout<<"proton: "<<event_p<<" neutron: "<<event_n<<" ratio: "<<ratio<<" stats p: "<<stats_p<<" stats n: "<<stats_n<<endl;



//code outline:
//use the SRC and MF cuts found previously

//only look at 4.4 GeV data for now

//correct both n and p for efficiencies (maps and parameterization
//calc 3He(e,e'n)/sigN / 3He(e,e'p)/sigP
//calc 4He(e,e'n)/sigN / 4He(e,e'p)/sigP --should be 1

//calc 12C(e,e'n)/sigN / 12C(e,e'p)/sigP --should be 1 (n det eff from 3He data)



}
