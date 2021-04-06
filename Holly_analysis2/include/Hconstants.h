#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__


const double Me      = .000511;         // Electron     mass in GeV
const double Mp      = 0.93827;         // Proton       mass in GeV
const double Mn      = 0.93957;         // Neutron      mass in GeV
const double Mpip    = 0.13957;         // Charged pion mass in GeV
const double He3_be  = 0.00770;         // 3He binding energy in GeV
const double M_tar   = 2*Mp+Mn-He3_be;  // 3Helium  mass in GeV (target)
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
const double amu     = 0.931494;
const double c_cm_ns = 29.9792458;
const double ns2s    = 1.0E-09 ;

//neutron background correction for neutron efficiency
const double corr_pp[3] = {0.731444,0.665635,0.318146};//bg subtraction for e,e'pp
const double corr_pppipi[2] = {0.75128,0.517747};//bg subtraction for e,e'pppipi


//SRC cuts:
const double vtx_min_C = 4.0;//cm
const double vtx_max_C = 7.0;//cm
const double vtx_min_He3 = -4.0;//cm
const double vtx_max_He3 = 0.5;//cm
const double vtx_min_He4 = -3.5;//cm
const double vtx_max_He4 = 2.5;//cm
const double xb_cut = 1.2;
const double thetapq_cut = 25.0;//deg
const double pq_min = 0.62;
const double pq_max_n = 1.1;
const double pq_max_p = 0.96;
const double mmiss_nomCut = 1.1;
const double pmiss_nomCut = 0.3;
const double Q2_min = 0.0;
//SRC neutron/smeared cuts:
const double mmiss_nomCut_n = 1.13;//1.1;//1.13;
const double pmissSRC_nomCut_n = 0.32;//0.3;//0.32;


//MF cuts:
const double y_cut_min = -0.05;
const double y_cut_max = 0.2;
const double nu_min = 0.9;
const double nu_max = 1.6;
const double thetapq_MFcut = 7.0;//deg
const double emiss_nomCut = 0.08;
const double pmiss_MFnomCut = 0.22;
const double emiss_nomCut_n = 0.25;//0.31;
const double pmiss_MFnomCut_n = 0.295;// 0.25;

#endif
