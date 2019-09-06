{
  gSystem->Load("libTELoss.so");

  // Creating instance of Energy Loss
  TCLASCorrectionsG11A  fCORR;
  //
  // The vector to be corrected 
  // and the vertex !
  //
  TVector3 v3_p(-0.8,-.25,0.2);
  TVector3 v3_v(0.0,0.0,-12.);

  TLorentzVector  vP;
  vP.SetVectM(v3_p,0.938);
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
  TLorentzVector vPc = fCORR.GetVector(vP,-1,2212);
  printf("x y z mag = %9.5f  %9.5f  %9.5f  %9.5f\n",
	 vP.Px(),vP.Py(),vP.Pz(),vP.P());
  printf("x y z mag = %9.5f  %9.5f  %9.5f  %9.5f\n",
	 vPc.Px(),vPc.Py(),vPc.Pz(),vPc.P());

}
