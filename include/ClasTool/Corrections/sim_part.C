{


  double   pmass[3] = {0.135,0.4973,0.938};
  TRandom  fRand;

  FILE *fp = fopen("evnt_data.evt","w");
  for(int i=0;i<20000;i++){
    double  z_targ = fRand.Rndm()*40.-10.;
    int     part   = (int) (fRand.Rndm()*3.);
    double  p_mom  = (fRand.Rndm()*1.2)+0.130;
    double  p_th   = (fRand.Rndm()*180.)/57.29;
    double  p_phi  = (fRand.Rndm()*360.)/57.29;
    TVector3 v3p;
    v3p.SetMagThetaPhi(p_mom,p_th,p_phi);
    fprintf(fp,"%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	   pmass[part],0.,0.,z_targ,
	   v3p.x(),v3p.y(),v3p.z());
  }
  fclose(fp);
}
