#include "TIdentificator.h"
#include "TClasTool.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"
#include "TMath.h"

int main(int argc, char **argv)
{
  TClasTool *input = new TClasTool();  
  input->InitDSTReader("ROOTDSTR");
    
  if(argc == 1) {
    char File[200];
    system("ls -1 *.root > dataFiles.txt");
    ifstream in("dataFiles.txt", ios::in);
    if (!in) {
      cerr << "File Not Opened!" << endl;
      exit(1);
    }
    while (in >> File) {
      input->Add(File);
    }
    in.close();
    system("rm dataFiles.txt");
  } 
  else 
    {
      for(int kklk=1; kklk<argc; kklk++)
	input->Add(argv[kklk]);
    }
  
  TIdentificator *t = new TIdentificator(input);
  TFile *output = new TFile("particle_data.root", "RECREATE", "Data of the Tree");
  TTree *tree = new TTree("data", "Tree that holds the data");

  //Reconstructed Variables
  //Write information for all events with at least 1 entry in EVNT Bank
  Int_t gPart;
  Int_t CCPart;
  Int_t DCPart;
  Int_t ECPart;
  Int_t SCPart;
  Float_t Stat[50];
  Int_t StatCC[50];
  Int_t StatDC[50];
  Int_t StatEC[50];
  Int_t StatSC[50];
  Int_t DC_Status[50];
  Int_t SC_Status[50];
  Int_t SC_Pad[50];

  Float_t EC_in[50];
  Float_t EC_out[50];
  Float_t EC_tot[50];
  Float_t EC_U[50];
  Float_t EC_V[50];
  Float_t EC_W[50];
  Float_t EC_X[50];
  Float_t EC_Y[50];
  Float_t EC_Z[50];
  Float_t Nphe[50];
  Float_t SC_Time[50];
  Float_t SC_Path[50];
  Float_t CC_Time[50];
  Float_t CC_Path[50];
  Float_t CC_Chi2[50];
  Float_t EC_Time[50];
  Float_t EC_Path[50];

  Float_t charge[50];
  Int_t id[50]; //Id in EVNT bank
  Float_t beta[50];
  Float_t mass[50];
  Float_t p[50];
  Float_t px[50];
  Float_t py[50];
  Float_t pz[50];
  Float_t theta[50];
  Float_t phi[50];
  Float_t z[50];
  Float_t theta_pq[50];

  Int_t NRun;
  Float_t STT;

  Float_t Xb;
  Float_t Q2;
  Float_t W;
  Float_t Nu;
  Float_t Yb;

  Int_t num_g;
  Int_t id_g[50];
  Float_t p_g[50];
  Float_t px_g[50];
  Float_t py_g[50];
  Float_t pz_g[50];
  Float_t theta_g[50];
  Float_t phi_g[50];
  Float_t z_g[50];

  // Load initial event
  input->Next();
  // This is crucial. Subsequent loading will occur in for loop increment

  // Determine if this is a simulation run based on the number of generated particles in the first event
  num_g = input->GetNRows("GSIM");
  const bool sim_run = (num_g>0);
  if (sim_run)
    {
      cout << "This is a simulation run. Generator branches will be saved.\n";
    }
  else
    {
      cout << "This is not a simulation run. No generator information.\n";
    }

  //Reconstructed Branches
  tree->Branch("gPart",&gPart,"gPart/I");
  tree->Branch("CCPart",&CCPart,"CCPart/I");
  tree->Branch("DCPart",&DCPart,"DCPart/I");
  tree->Branch("ECPart",&ECPart,"ECPart/I");
  tree->Branch("SCPart",&SCPart,"SCPart/I");
  tree->Branch("Stat",Stat,"Stat[gPart]/F");
  tree->Branch("StatCC",StatCC,"StatCC[gPart]/I");
  tree->Branch("StatDC",StatDC,"StatDC[gPart]/I");
  tree->Branch("StatEC",StatEC,"StatEC[gPart]/I");
  tree->Branch("StatSC",StatSC,"StatSC[gPart]/I");
  tree->Branch("DC_Status",DC_Status,"DC_Status[gPart]/I");
  tree->Branch("SC_Status",SC_Status,"SC_Status[gPart]/I");
  tree->Branch("SC_Pad",SC_Pad,"SC_Pad[gPart]/I");
  tree->Branch("EC_in",EC_in,"EC_in[gPart]/F");
  tree->Branch("EC_out",EC_out,"EC_out[gPart]/F");
  tree->Branch("EC_tot",EC_tot,"EC_tot[gPart]/F");
  tree->Branch("EC_U",EC_U,"EC_U[gPart]/F");
  tree->Branch("EC_V",EC_V,"EC_V[gPart]/F");
  tree->Branch("EC_W",EC_W,"EC_W[gPart]/F");
  tree->Branch("EC_X",EC_X,"EC_X[gPart]/F");
  tree->Branch("EC_Y",EC_Y,"EC_Y[gPart]/F");
  tree->Branch("EC_Z",EC_Z,"EC_Z[gPart]/F");
  tree->Branch("Nphe",Nphe,"Nphe[gPart]/F");
  tree->Branch("SC_Time",SC_Time,"SC_Time[gPart]/F");
  tree->Branch("SC_Path",SC_Path,"SC_Path[gPart]/F");
  tree->Branch("EC_Time",EC_Time,"EC_Time[gPart]/F");
  tree->Branch("EC_Path",EC_Path,"EC_Path[gPart]/F");
  tree->Branch("CC_Time",CC_Time,"CC_Time[gPart]/F");
  tree->Branch("CC_Path",CC_Path,"CC_Path[gPart]/F");
  tree->Branch("CC_Chi2",CC_Chi2,"CC_Chi2[gPart]/F");
  tree->Branch("Charge",charge,"charge[gPart]/F");
  tree->Branch("particle",id, "id[gPart]/I");
  tree->Branch("Beta",beta,"beta[gPart]/F");
  tree->Branch("Mass",mass,"mass[gPart]/F");
  tree->Branch("Momentum",p,"p[gPart]/F");
  tree->Branch("Momentumx",px,"px[gPart]/F");
  tree->Branch("Momentumy",py,"py[gPart]/F");
  tree->Branch("Momentumz",pz,"pz[gPart]/F");
  tree->Branch("Theta",theta,"theta[gPart]/F");
  tree->Branch("Phi",phi,"phi[gPart]/F");
  tree->Branch("TargetZ",z,"z[gPart]/F");
  tree->Branch("Thetapq",theta_pq,"theta_pq[gPart]/F");
  tree->Branch("NRun",&NRun,"NRun/I");
  tree->Branch("STT",&STT,"STT/F");
  tree->Branch("BjorkenX",&Xb,"Xb/F");
  tree->Branch("Q2",&Q2,"Q2/F");
  tree->Branch("W",&W,"W/F");
  tree->Branch("Nu",&Nu,"Nu/F");
  tree->Branch("Yb",&Yb,"Yb/F");

  if (sim_run)
    {
      tree->Branch("Number_g",&num_g,"num_g/I");
      tree->Branch("particle_g",&id_g,"id_g[num_g]/I");
      tree->Branch("Momentum_g",p_g,"p_g[num_g]/F");
      tree->Branch("Momentumx_g",px_g,"px_g[num_g]/F");
      tree->Branch("Momentumy_g",py_g,"py_g[num_g]/F");
      tree->Branch("Momentumz_g",pz_g,"pz_g[num_g]/F");
      tree->Branch("Theta_g",theta_g,"theta_g[num_g]/F");
      tree->Branch("Phi_g",phi_g,"phi_g[num_g]/F");
      tree->Branch("TargetZ_g",z_g,"z_g[num_g]/F");
    } 

  Int_t nEntries = input->GetEntries();
  cout << "Total Entries = "<<nEntries<<endl;

  // Loop over events, loading each through Next() command
  for (int event = 0; event < nEntries; event++, input->Next()) {
    
    if(event%10000 == 0) cout << "Events Processed: "<< event << endl;

    //Reconstucted Variables
    gPart = t->gPart();
    CCPart = t->CCPart();
    DCPart = t->DCPart();
    ECPart = t->ECPart();
    SCPart = t->SCPart();

    // Fill once per event
    if (gPart > 0)
      {
	NRun = t->NRun();
	STT = t->STT();
	Xb = t->Xb();
	Q2 = t->Q2();
	W = t->W();
	Nu = t->Nu();
	Yb = t->Yb();
      }

    // Loop over reconstructed particles
    for (Int_t i = 0; i < TMath::Min(gPart, 50); i++)
      {
	Stat[i] = t->Status(i);
	StatCC[i] = t->StatCC(i);
	StatDC[i] = t->StatDC(i);
	StatEC[i] = t->StatEC(i);
	StatSC[i] = t->StatSC(i);
	DC_Status[i] = t->DCStatus(i);
	SC_Status[i] = t->SCStatus(i);
	SC_Pad[i] = t->SCpdht(i);
	
	EC_in[i] = t->Ein(i);
	EC_out[i] = t->Eout(i);
	EC_tot[i] = t->Etot(i);
	EC_U[i] = t->EC_U(i);
	EC_V[i] = t->EC_V(i);
	EC_W[i] = t->EC_W(i);
	EC_X[i] = t->EC_X(i);
	EC_Y[i] = t->EC_Y(i);
	EC_Z[i] = t->EC_Z(i);
	Nphe[i] = t->Nphe(i);
	SC_Time[i] = t->TimeSC(i);
	SC_Path[i] = t->PathSC(i);
	CC_Chi2[i] = t->Chi2CC(i);
	CC_Time[i] = t->TimeCC(i);
	CC_Path[i] = t->PathCC(i);
	EC_Time[i] = t->EC_time(i);
	EC_Path[i] = t->EC_path(i);
	charge[i] = t->Charge(i);
	id[i] = t->Id(i);
	beta[i] = t->Betta(i);
	mass[i] = TMath::Sqrt(t->Mass2(i));
	p[i] = t->Momentum(i);
	px[i] = t->Px(i);
	py[i] = t->Py(i);
	pz[i] = t->Pz(i);
	theta[i]= t->FidTheta(i);
	phi[i]= t->FidPhi(i);
	z[i] = t->Z(i);
	theta_pq[i] = (TMath::ACos(t->CosThetaPQ(i)) * TMath::RadToDeg());
      }
    
    // Loop over simulated particles
    num_g = input->GetNRows("GSIM");
    for (int j=0 ; j<num_g ; j++)
      {
	id_g[j]=t->Id(j,1);
	p_g[j]=t->Momentum(j,1);
	px_g[j]=t->Px(j,1);
	py_g[j]=t->Py(j,1);
	pz_g[j]=t->Pz(j,1);
	theta_g[j] = t->FidTheta(j,1);
	phi_g[j] = t->FidPhi(j,1);
	z_g[j] = t->Z(j,1);
      }

    if ((sim_run && (num_g > 0)) || ((!sim_run) && (gPart>0) && (gPart<50)))
      tree->Fill();
  }
  
  output->Write();
  output->Close();
  cout << "Done." << endl;
  
  return 0;
}
