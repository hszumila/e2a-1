# Install the calibration data so that programs can read it.
INSTALL (FILES
	el_Epratio_mom_2261.root el_Epratio_mom_4461.root el_Epratio_mom_1161.root
	el_Epratio_mom_4461_sim.root
	FCP_2261_2250.dat FCP_4461_2250.dat FCP_1161_1500.dat FCP_1161_750.dat
	EMCP_2261_2250.dat EMCP_4461_2250.dat EMCP_1161_1500.dat EMCP_1161_750.dat
	protdeltat_mom_2261_2250.root protdeltat_mom_4461_2250.root protdeltat_mom_1161_1500.root protdeltat_mom_1161_750.root
	pipdeltat_mom_2261_2250.root pipdeltat_mom_4461_2250.root pipdeltat_mom_1161_1500.root pipdeltat_mom_1161_750.root
	pimdeltat_mom_2261_2250.root pimdeltat_mom_4461_2250.root pimdeltat_mom_1161_1500.root pimdeltat_mom_1161_750.root
	PFID_2261_2250.dat PFID_4461_2250.dat PFID_1161_1500.dat PFID_1161_750.dat
	PIMFID_1161_1500.dat PIMFID_1161_750.dat
	run_table.dat
	prot_mom_corr_18338_18438.dat
	vz_4461_4He.dat vz_4461_12C.dat vz_4461_3He.dat vz_1161_3He.dat
	vz_2261_4He.dat vz_2261_12C.dat vz_2261_3He.dat vz_2261_56Fe.dat
	n_pathlength_corr_4461.dat
        DESTINATION "$ENV{E2A_INSTALL}")


#	These were the old root file vertex corrections. Now replaced with 2-param text files.
#	vz_4461_4He.root vz_4461_12C.root vz_1161_3He.root
#	vz_2261_4He.root vz_2261_12C.root vz_2261_3He.root vz_2261_56Fe.root
