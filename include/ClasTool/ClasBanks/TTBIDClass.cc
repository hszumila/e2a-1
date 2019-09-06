#include "TTBIDClass.h"
ClassImp(TTBIDClass)

TTBIDClass::TTBIDClass(TTBIDClass *TmpTBID)
{
	track = TmpTBID->track;
	sec   = TmpTBID->sec;
	beta  = TmpTBID->beta;
	vtime = TmpTBID->vtime;
	
	sc_stat = TmpTBID->sc_stat;
	sc_id   = TmpTBID->sc_id;
	sc_time = TmpTBID->sc_time;
	sc_qual = TmpTBID->sc_qual;
	sc_vtime = TmpTBID->sc_vtime;
	sc_beta  = TmpTBID->sc_beta;

	cc_stat  = TmpTBID->cc_stat;
	cc_id    = TmpTBID->cc_id;
	cc_time  = TmpTBID->cc_time;
	cc_qual  = TmpTBID->cc_qual;
	cc_vtime = TmpTBID->cc_vtime;
	cc_beta  = TmpTBID->cc_beta;

	ec_stat  = TmpTBID->ec_stat;
	ec_id    = TmpTBID->ec_id;
	ec_time  = TmpTBID->ec_time;
	ec_qual  = TmpTBID->ec_qual;
	ec_vtime = TmpTBID->ec_vtime;
	ec_beta  = TmpTBID->ec_beta;

	st_stat  = TmpTBID->st_stat;
	st_id    = TmpTBID->st_id;
	st_time  = TmpTBID->st_time;
	st_qual  = TmpTBID->st_qual;
	st_vtime = TmpTBID->st_vtime;
	st_beta  = TmpTBID->st_beta;

	lac_stat = TmpTBID->lac_stat;
	lac_id   = TmpTBID->lac_id;
	lac_time = TmpTBID->lac_time;
	lac_qual = TmpTBID->lac_qual;
	lac_vtime = TmpTBID->lac_vtime;
	lac_beta = TmpTBID->lac_beta;
}

void TTBIDClass::Print(void)
{

}
/*
void TTBIDClass::Print(void)
{
  cout << "TTBERClass ::" << endl;
  cout << "q_over_p   = " << q_over_p << endl;
  cout << "lambda     = " << lambda  << endl;
  cout << "phi        = " << phi     << endl;
  cout << "d0         = " << d0      << endl;
  cout << "z0         = " << z0      << endl;
  cout << "chi2       = " << chi2    << endl;
  cout << "layinfo1   = " << layinfo1<< endl;
  cout << "layinfo2   = " << layinfo2<< endl;
  cout << "Error Matrix:" << endl;
}*/
