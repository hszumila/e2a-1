
TString TIdentificator::GetCategorization(Int_t k)
{
  Int_t number_dc = fCT->GetNRows("DCPB");
  Int_t number_cc = fCT->GetNRows("CCPB");
  Int_t number_sc = fCT->GetNRows("SCPB");
  Int_t number_ec = fCT->GetNRows("ECPB");
  
  TString partId;
  
  partId = "not recognized";
  
  if (number_dc != 0) {
    if (k == 0 &&
	Status(0) > 0 && Status(0) < 100 &&
	Charge(0) == -1 &&
	number_cc != 0 && number_ec != 0 && number_sc != 0 &&
	StatCC(0) > 0 && StatSC(0) > 0 &&
	StatDC(0) > 0 && StatEC(0) > 0 &&
	DCStatus(0) > 0 &&
	Nphe(0) > 25 &&
	
	// EC Cuts Used By The CT Group (approved by the CLAS NPWG)
	EC_U(0) > 40 &&
	EC_V(0) < 360 &&
	EC_W(0) < 395 && 
	Ein(0)   > 0.05 &&
	Eout(0)  > 0.01 &&
	Eout(0)/Momentum(0) > ( 0.271*(1-0.3/sqrt( 0.5                      )) - Ein(0)/Momentum(0) )  &&
	Eout(0)/Momentum(0) < ( 0.271*(1+0.3/sqrt( electronPID_b(Momentum(0)) )) - Ein(0)/Momentum(0) ) ) {
      partId = "electronNoFid";
      if (FidCheckCut() == 1)
	partId = "electron";
    }
    else if (k == 0 &&
	     Status(0) > 0 && Status(0) < 100 &&
	     Charge(0) == -1 &&
	     number_cc != 0 && number_ec != 0 && number_sc != 0 &&
	     StatCC(0) > 0 && StatSC(0) > 0 &&
	     StatDC(0) > 0 && StatEC(0) > 0 &&
	     DCStatus(0) > 0){
      partId = "electronNoPid";
    }
    

      if (k > 0) {
	
	//Photons
	if (Charge(k)==0 && Betta(k)>0.95 && Betta(k)<1.95 && StatEC(k)>0)
	  partId = "photon";
	//Neutrons
	if (Charge(k)==0 && Betta(k)<0.95 && StatEC(k)>0)
          partId = "neutron";
	
	//Positive Particles
	if (Charge(k) == 1 &&
	    Status(k) > 0 && Status(k) < 100 &&
	    StatDC(k) > 0 && DCStatus(k) > 0) {
	  partId = "positive";

	  //Pions
	  if (Momentum(k)>=2.7 &&
	      number_cc != 0 && StatCC(k) > 0 &&
	      Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3){
	    partId = "piplusNoFid";
	    if (FidCheckCutPiPlus(k) == 1)
	      partId = "piplus";
	  }
	  else if( Momentum(k) < 2.7 && 
		   number_sc != 0 && StatSC(k) > 0 &&
		   ((Momentum(k) < 1 &&
		     TimeCorr4(0.139,k) >= -1.46 &&
		     TimeCorr4(0.139,k) <= 0.15) ||
		    (Momentum(k) >=1 &&
		     TimeCorr4(0.139,k) >= -1.38 &&
		     TimeCorr4(0.139,k) <= 0.53)) ){
	    partId ="piplusNoFid";
	    if (FidCheckCutPiPlus(k) == 1)
              partId = "piplus";
	  }
	  
	  //Protons
	  if (Momentum(k) < 2.8 && 
	      number_sc != 0 && StatSC(k) > 0 &&
	      StatSC(0) > 0 && SCStatus(0) > 29 && SCStatus(k) > 29 && //Added by Barak (August 2015)
	      ProtonPID(TimeCorr4(0.938,k), Momentum(k))){
	    partId = "protonNoFid";
	    if (FidCheckCutPiPlus(k) == 1)
	      partId = "proton";
	  }
	
	  //Positrons
	  if (number_cc != 0 && number_ec != 0 && number_sc != 0 &&
	      StatCC(k) > 0 && StatSC(k) > 0 && StatEC(k) > 0 &&
	      Nphe(k) > 25 &&
	      Etot(k) / 0.27 + 0.4 > Momentum(k) &&
	      Etot(k) / 0.27 - 0.4 < Momentum(k) &&
	      Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
	      Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k)){
	    partId = "positron";
	  }
	}
	
	//Negative Pions
	if (Charge(k) == -1 &&
	    Status(k) > 0 && Status(k) < 100 &&
	    StatDC(k) > 0 && DCStatus(k) > 0 && Momentum(k)>=2.7 &&
	    number_cc != 0 && StatCC(k) > 0 &&
	    Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3){
	  partId = "piminusNoFid";
	  if (FidCheckCutPiMinus(k) == 1)
	    partId = "piminus";
	}
	else if (Charge(k) == -1 &&
		 Status(k) > 0 && Status(k) < 100 &&
		 StatDC(k) > 0 && DCStatus(k) > 0 &&
		 Momentum(k) < 2.7 &&
		 number_sc != 0 && StatSC(k) > 0 &&
		 ((Momentum(k) < 1 &&
		   TimeCorr4(0.139,k) >= -1.46 &&
		   TimeCorr4(0.139,k) <= 0.15) ||
		  (Momentum(k) >=1 &&
		   TimeCorr4(0.139,k) >= -1.38 &&
		   TimeCorr4(0.139,k) <= 0.53)) ){
	  partId = "piminusNoFid";
	  if (FidCheckCutPiMinus(k) == 1)
            partId = "piminus";
	} 
	
      }
  }
  
  return partId;
}


TString* TIdentificator::GetCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");

    if (fPartIds != 0) delete [] fPartIds;
    fPartIds = new TString[number];

    if (number != 0) {
        for (Int_t i = 0; i < number; i++)
            fPartIds[i] = GetCategorization(i);
    }

    return fPartIds;
}



void TIdentificator::PrintCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");
    fPartIds = GetCategorization();

    if (fPartIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(fPartIds+i) << endl;
        cout << endl;
    }
}



void TIdentificator::PrintCategorization(TString* partIds)
{
    Int_t number = fCT->GetNRows("EVNT");

    if (partIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(partIds+i) << endl;
        cout << endl;
    }
}

bool TIdentificator::ProtonPID(Double_t time, Double_t p){
  
  Double_t up1[10] = {120.251, -1168.19, 5237.42, -13638.1, 22325.3, -23649.6, 16177.8, -6894.16, 1663.4, -173.474};
  Double_t up2[10] = {1.59223, -1.49056, 0.676338, -0.104644, 0, 0, 0, 0, 0, 0};
  Double_t down1[10] = {-26.8257, 153.155, -408.979, 673.476, -790.433, 708.255, -474.843, 217.295, -58.6269, 6.90981};
  Double_t down2[10] = {-1.1009, 0.719988, -0.280016, 0.0319352, 0, 0, 0, 0, 0, 0};
  
  Double_t time_up=0;
  Double_t time_down=0;
  
  for(int i=0; i<10; i++){
    time_up += (p<0.8)?(up1[i]*pow(p,i)):(up2[i]*pow(p,i));
    time_down += (p<0.8)?(down1[i]*pow(p,i)):(down2[i]*pow(p,i));
  }
  
  return (time>time_down && time<time_up);
}

Double_t TIdentificator::electronPID_b(Double_t p){
  
  if (p<0.7) return 0.85;
  
  Double_t b[22] = {0.85, 0.8, 0.85, 1.05, 1.1, 1.35, 1.35, 1.45, 1.35, 1.35, 1.35, 1.3, 1.35, 1.35, 1.5, 1.6, 1.8, 1.8, 1.8, 1.8,1.8,1.8};
  return b[ ((int)( (p-0.5)/0.2 )) ];
}

