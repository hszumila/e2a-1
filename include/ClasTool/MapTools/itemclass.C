{

gROOT->Reset();

TSQLServer *db = TSQLServer::Connect("mysql://quark","root","");
db->SelectDataBase("calib"); 

TItemClass T;
T->SetItem(db,"RUN_CONTROL","beam","energy");
if(T->FindItemEntry(db,14567,"suleiman")!=0){
	cout << "ERROR" << endl;
}

}