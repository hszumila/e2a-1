#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

double max_outta_three ( double A , double B, double C ){
  
  double temp_max, max_val;
  
  temp_max = TMath::Max(A,B       );
  max_val  = TMath::Max(C,temp_max);
  
  return max_val;
}

TLorentzVector correctPn(TLorentzVector mom_n, double eIn, double eOut){
  double a,b,c;
  
  TVector3 mom = mom_n.Vect();
  double P = mom.Mag();
  
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    a = 0.230042;
    b = -0.3888288;
    c = 0.1421874;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    a = 0.5843257;
    b = -0.8823494;
    c = 0.309747;
   }
  //both
  else {
    a = 0.230042;
    b = -0.3888288;
    c = 0.1421874;
  }

  double newMag;
  if (P>2.1){newMag = P - P*(a+b*2.1+c*2.1*2.1);}
  else if (P<0.5){newMag = P - P*(a+b*0.5+c*0.5*0.5);}
  else{newMag = P - P*(a+b*P+c*P*P);}
  
  TVector3 unit = mom.Unit();
  unit.SetMag(newMag);

  mom_n.SetVect(unit);
  
  return mom_n;

  
}

TVector3 correctPnVec(TVector3 mom, double eIn, double eOut){
  double a,b,c;
  
  double P = mom.Mag();
  
  //inner EC
  if (eIn>0.01 && eOut<0.01){
    a = 0.230042;
    b = -0.3888288;
    c = 0.1421874;
  }
  //outer EC
  else if (eIn<0.01 && eOut>0.01){
    a = 0.5843257;
    b = -0.8823494;
    c = 0.309747;
   }
  //both
  else {
    a = 0.230042;
    b = -0.3888288;
    c = 0.1421874;
  }

  double newMag;
  if (P>2.1){newMag = P - P*(a+b*2.1+c*2.1*2.1);}
  else if (P<0.5){newMag = P - P*(a+b*0.5+c*0.5*0.5);}
  else{newMag = P - P*(a+b*P+c*P*P);}
  
  TVector3 unit = mom.Unit();
  unit.SetMag(newMag);
  
  return unit;
  
}

TLorentzVector calculatePn(double beta, TLorentzVector mom){

  TVector3 pn = mom.Vect();
  TVector3 unit = pn.Unit();
  double corrP = Mn*beta/sqrt(1.0-beta*beta);
  unit.SetMag(corrP);
  mom.SetVect(unit);
  return mom;

}


bool p_points_to_EC_fid(TVector3 pm, double vtx_z, double dist){
	double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;
	double pi = 22./7.;

	const double ec_theta =  65.*pi/180.; // degrees
	const double z0       =  516.9152;//563.1; // cm

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


Double_t func_fit(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t value;
  if     (xx <  par[0]){value = par[1]*xx    +par[2];}
  else if(xx >= par[0]){value = par[1]*par[0]+par[2];}
  return value;
}

Double_t func_fitEG2(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t value;
  if     (xx <  par[0]){value = par[1]*xx    +par[2];}
  else if(xx >= par[0]){value = par[1]*par[0]+par[2];}


  if(xx>0.5 && xx<1.){value *= 1.24;}
  else if (xx>=1. && xx<1.5) {value *= 1.21;}
  else if (xx>=1.5 && xx<2.) {value *= 1.22;}
  else if (xx>=2. && xx<2.5) {value *= 1.15;}
  
  return value;
}

Double_t func_fitExt(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t value;
  value = par[0]*xx+par[1];
  //if     (xx <  par[0]){value = par[1]*xx    +par[2];}
  //else { value = par[3]*xx+(par[1]-par[3])*par[0]+par[2];}
  //else if(xx >= par[0] && xx<par[3]){value = par[4]*xx+(par[1]-par[4])*par[0]+par[2];}
  //else {value = par[4]*par[3]+(par[1]-par[4])*par[0]+par[2];}
  return value;
}
// ====================================================================================================================================================
void fit_gfunction(TH1F * gPlot, float minx, float maxx){
  Double_t par[3] = {1.8,0.3,-0.1};
  TF1* Func = new TF1("Func", func_fit, 0.6,2.,3);
  Func -> SetParameters(par);
  Func -> SetParLimits(0,1.55,1.8);
  gPlot -> Fit("Func","EM","", minx, maxx);
}

void fit_gfunctionExt(TH1F * gPlot, float minx, float maxx){
  //Double_t par[5] = {1.1,0.2,0.1,1.6,0.35};
  Double_t par[2] = {1.1,0.2};
  TF1* Func = new TF1("Func", func_fitExt, 0.7,2.1,2);
  Func -> SetParameters(par);
  //Func -> SetParLimits(0,1.1,1.5);
  //Func -> SetParLimits(3,1.65,1.7);
  gPlot -> Fit("Func","EM","", minx, maxx);
}

void fit_gfunctionEG2(TH1F * gPlot, float minx, float maxx){
  Double_t par[3] = {1.57,0.28,-0.11};
  TF1* Func = new TF1("Func", func_fitEG2, 0.6,2.1,3);
  //Func -> SetParameters(par);
  Func->FixParameter(0,1.57);
  Func->FixParameter(1,0.28);
  Func->FixParameter(2,-0.11);
  //Func -> SetParLimits(0,0.9,1.3);
  //Func -> SetParLimits(3,1.55,1.65);
  gPlot -> Fit("Func","EM","", minx, maxx);
}
// ====================================================================================================================================================
// polynomial background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}
// Lorentzian Peak function
Double_t gaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*(pow((x[0]-par[1])/par[2],2.0)));  
}
// Sum of background and peak function
Double_t fitSB(Double_t *x, Double_t *par) {
  return background(x,par) + gaussianPeak(x,&par[4]);
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
  
  //4.4 GeV 3He:
  double ec_extr_res = 0.106332;
  double dt = 0.387394;          //ec_time_res;
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
  //my best calc (Oct 2020):
  /*
  if (pm<1.65){value = 0.415133*pm -0.188579;}
  else {value = 0.415133*1.65 -0.188579;}
  return value;
  */
  //my new calc with plc (March 2021):
 // if (pm<1.7){value = 0.368212*pm -0.154081;}
 // else if(pm>=1.7) {value = 0.368212*1.7 -0.154081;}//0.37129*1.65 -0.156182;}
  //return value;
  //if (pm<1.65){value = 0.95*(0.371291*pm -0.156184);}
    //else if(pm>=1.65) {value = 0.95*(0.371291*1.65 -0.156184);}//0.37129*1.65 -0.156182;}
   // return value;
  //if (pm<1.8){value = 0.363075*pm -0.150405;}
  //else if(pm>=1.8) {value = 0.363075*1.8 -0.150405;}//0.37129*1.65 -0.156182;}
  //return value;
  //double params[5] = {0.932139, 0.176687, -0.236691, 1.7, 0.549818};
  //double params[5] = {1.3, 0.289927, -0.0986945, 1.65, 0.858808};
  //if     (pm < params[0]){value = params[1]*pm + params[2];}
  //else if(pm>= params[0] && pm<params[3]){value = params[4]*pm+(params[1]-params[4])*params[0] + params[2];}
  //else {value = params[4]*params[3]+(params[1]-params[4])*params[0]+params[2];}
  //return value;  
    
  //Meytal's:
  if (pm<1.57&&pm>0.5){value = 0.28*pm -0.11;}
  else if (pm<=0.5){value = 0.28*0.5 -0.11;}
  else {value = 0.28*1.57 -0.11;}
    
  if (0.5<=pm && pm<1){value = 1.24*value;}
  else if (1<=pm && pm<1.5){value = 1.21*value;}
  else if (1.5<=pm && pm<2){value = 1.22*value;}
  else{value = 1.15*value;}
  return value;
  
  
}


#endif
