#include "TTBTRClass.h"
ClassImp(TTBTRClass)

TTBTRClass::TTBTRClass(TTBTRClass *TmpTBTR)
{
	x = TmpTBTR->x;
	y   = TmpTBTR->y;
	z  = TmpTBTR->z;
	px = TmpTBTR->px;
	py = TmpTBTR->py;
	pz = TmpTBTR->pz;
	q = TmpTBTR->q;
	chi2 = TmpTBTR->chi2;
	itr_sec = TmpTBTR->itr_sec;
	itr_hbt = TmpTBTR->itr_hbt;
}

void TTBTRClass::Print(void)
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
