//------------------------------------------------------------------ 
//
//  This example shows how to get Run Energy from 
//  Calibration Database. (mySQL).
//  To make it work on your system (if you're not runnning this code
//  at Jefferson Lab) just change hostname and username variables
//
//-------------------------------------------------------------------
//  Author : Gagik Gavalian (UNH), Date : 06/10/2000 
//-------------------------------------------------------------------

{

gROOT->Reset();

//TMapUtils  CalibMap("alquds.jlab.org","calib","dbmanager");
TMapUtils  CalibMap("quark","calibnew","dbuser");
CalibMap->PrintInfo();
CalibMap->SetDebugMode(0x20);
TArrayF  BeamEnergy;

CalibMap->Get_Map_Float("RUN_CONTROL","beam","energy",17870,&BeamEnergy);

CalibMap->PrintFloatArray(BeamEnergy);

cout << " Beam Energy = "  << endl;
}
