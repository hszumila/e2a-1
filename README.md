# e2a
Repository for everything e2a. In this repository, you will find:

    - skim_tree.cpp : this is the flagship code of the repository
                         it skims data for a selected reaction type, applying PID and 
			 fiducial cuts
    - cut_library : a library with classes that can apply PID, Fiducial cuts,
                         corrections, etc.
    - calibration_data : a directory with data files read in by the cut_library
    - cut_production : has programs used to generate the data files in calibration_data
    - maps : a directory containing acceptance maps, and the code used to create them
    - write_tree : a directory for code used to convert BOS files into a root format
                      suitable for skimming
    - Rey_neutrons : a directory for Rey's e2a work
    - Axel_analysis : a directory for Axel's e2a work
    - Adin_analysis : a directory for Adin's e2a work

For a minimal working example of how one might analyze skimmed data,
    see the program Axel_analysis/example_analysis.cpp

Features we would like to add in the future:
	 - Acceptance maps
	 - Resolution maps
	 - Background subtraction codes


URGENT THINGS THAT ARE MISSING
       - Vertex corrections for 4 GeV Iron
       - 1 GeV anything.