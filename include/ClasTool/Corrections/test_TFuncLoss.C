{
  gSystem->Load("libTELoss.so");

  TFuncLoss  fELoss;
  
  double pout[3],pin[3];

  pout[0] = 0.14;
  pout[1] = 0.15;
  pout[2] = 0.13;

  //  fELoss.Print(0);
  TVector3  aa;
  //  cout << aa << endl;
  fELoss.Loss(0,3.,0.938,pout,pin);


  
  printf("output momentum = %9.5f %9.5f %9.5f\n",pin[0],pin[1],pin[2]);
  
  for(double kk=0.05;kk<50.;kk+=0.05){
    //  printf("%9.5f  %9.5f\n",kk,fELoss.GetFunc(0,kk));
  }
  

}
