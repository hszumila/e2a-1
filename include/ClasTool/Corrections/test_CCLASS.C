{
  gSystem->Load("libTELoss.so");

  // Creating instance of Energy Loss
  TEnergyLoss  fEnL;
  //
  // This prints a warning and usage for the class
  //
  fEnL.PrintWarning();
  //
  // Initialize the Energy Loss for g11a data 
  //
  fEnL.Init(6);
  //
  // Set debug mode to print step by step corrections
  // To remove the printout use fEnL.SetDebugMode(0);
  //
  fEnL.SetDebugMode(0);
  FILE *fp = fopen("evnt_data.evt","r");
  //
  // The vector to be corrected 
  // and the vertex !
  //

  double v3p[3],v3v[3],mass;
  for(int i=0;i<20000;i++){

    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,
	   &v3v[0],&v3v[1],&v3v[2],
	   &v3p[0],&v3p[1],&v3p[2]);

  TVector3 v3_p(v3p[0],v3p[1],v3p[2]);
  TVector3 v3_v(v3v[0],v3v[1],v3v[2]);

  TVector3 v3_corr = fEnL.GetVector(v3_p,mass,v3_v,1,6);
  printf("%9.5f %9.5f %9.5f\n",mass,v3_p.Mag(),v3_corr.Mag());
  }
  fclose(fp);
}
