#include "TFuncLoss.h"


TFuncLoss::TFuncLoss(){

  fR_Spline[0]     = new TSpline3("SRLH2_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_RLH2,190);
  fR_Spline_REV[0] = new TSpline3("SRLH2_R",TEnergyLossNameSpace::DATA_RLH2,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[1]     = new TSpline3("SRSCINT_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_SCINT,190);
  fR_Spline_REV[1] = new TSpline3("SRSCINT_R",TEnergyLossNameSpace::DATA_SCINT,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[2]     = new TSpline3("SRCARBON_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_CARBON,190);
  fR_Spline_REV[2] = new TSpline3("SRCARBON_R",TEnergyLossNameSpace::DATA_CARBON,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[3]     = new TSpline3("SRLD2_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_LD2,190);
  fR_Spline_REV[3] = new TSpline3("SRLD2_R",TEnergyLossNameSpace::DATA_LD2,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[4]     = new TSpline3("SRLH3_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_RLH3,190);
  fR_Spline_REV[4] = new TSpline3("SRLH3_D",TEnergyLossNameSpace::DATA_RLH3,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[5]     = new TSpline3("SRAL_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_RAL,190);
  fR_Spline_REV[5] = new TSpline3("SRAL_R",TEnergyLossNameSpace::DATA_RAL,TEnergyLossNameSpace::DATA_POM,190);
  
  fR_Spline[6]     = new TSpline3("SRL4HE_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_L4HE,190);
  fR_Spline_REV[6] = new TSpline3("SRL4HE_R",TEnergyLossNameSpace::DATA_L4HE,TEnergyLossNameSpace::DATA_POM,190);

  fR_Spline[7]     = new TSpline3("SRAIR_D",TEnergyLossNameSpace::DATA_POM,TEnergyLossNameSpace::DATA_RAIR,190);
  fR_Spline_REV[7] = new TSpline3("SRAIR_R",TEnergyLossNameSpace::DATA_RAIR,TEnergyLossNameSpace::DATA_POM,190);

}

TFuncLoss::~TFuncLoss(){

}

int TFuncLoss::Loss(int material, double thikness, double pmass, double pout[], double pin[]){
  
  double    *ptr_Material = NULL;
  double  _part_mom = sqrt(pout[0]*pout[0]+pout[1]*pout[1]+pout[2]*pout[2]);
  double  poverm = _part_mom/pmass;

  pin[0] = pout[0];
  pin[1] = pout[1];
  pin[2] = pout[2];

  if(thikness<0.) return 1;
  if(material<0||material>7) return 0;
  if(poverm<0.05||poverm>50.) return -1;

  //  if(material==0) ptr_Material = DATA_RLH2;
  double  r_val = fR_Spline[material]->Eval(poverm);
  double  x  = r_val + thikness/pmass;
  double  piom = fR_Spline_REV[material]->Eval(x);

  double  pim  = piom*pmass;
  double  dp   = pim/_part_mom;

  pin[0] = pout[0]*dp;
  pin[1] = pout[1]*dp;
  pin[2] = pout[2]*dp;

  //  printf("%9.5f\n",r_val);
  return 1;
}

double  TFuncLoss::GetFunc(int material, double x_func)
{
  return  fR_Spline[material]->Eval(x_func);
}

void  TFuncLoss::Print(int material)
{
  double *ptr_todata=NULL;
  if(material==0){
    ptr_todata = TEnergyLossNameSpace::DATA_RLH2;
  }


  if(ptr_todata==NULL) return;
  for(int i=0;i<MATRIX_LEN;i++){
    printf("x : %9.5f %9.5f\n",TEnergyLossNameSpace::DATA_POM[i],ptr_todata[i]);
  }
}
