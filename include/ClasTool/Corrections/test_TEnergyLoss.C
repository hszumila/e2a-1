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
  fEnL.SetDebugMode(1);

  //
  // The vector to be corrected 
  // and the vertex !
  //
  TVector3 v3_p(.24,.25,0.23);
  TVector3 v3_v(0.0,0.0,-12.);
  //
  // Calling correction Routine
  // arguments :
  // 1 - particle vector
  // 2 - particle mass
  // 3 - particle vertex vector
  // 4 - the target material
  //     0 - empty target
  //     1 - LH2
  //     2 - LD2
  //     3 - L3He
  //     4 - L4He
  // 6 - target Cell Described in the Warning print out
  //
  //  Output: The corrected 3vector
  //
  TVector3 v3_corr = fEnL.GetVector(v3_p,0.938,v3_v,1,6);

}
