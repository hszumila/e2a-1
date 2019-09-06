{
gROOT->Reset();

TString slibpath(gSystem->Getenv("TOP_DIR"));
TString osname(gSystem->Getenv("OSCLAS"));

slibpath += "/slib/";
slibpath += osname.Data();
slibpath += "/";
TString MapLib   = slibpath + "libMapTools.so";
TString MySQLLib = slibpath + "libMySQL.so";

cout << "Loading Library : " << MapLib.Data() << endl;
gSystem->Load(MapLib.Data());
cout << "Loading Library : " << MySQLLib.Data() << endl;
gSystem->Load(MySQLLib.Data());
}
