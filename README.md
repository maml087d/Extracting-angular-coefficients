The purpose of this repository is to extract the angular coefficients from a EVNT file and generate the templates. The workwflow can be guided with the bash scripts in the main folder. This will automate the process of running the rivet analysis and then extracting the the relevant histograms from the .yoda file using the read_/*.py scripts in the rivet_anas folder.
Relevant bash scripts:
- extract_angcoeffs.sh
- maketemps_binned.sh
- setup_rivet.sh -> sets up a screen with the ATLAS environment

All of these use the ATLAS env: asetup 23.6.26,AthGeneration and automatically creates these environments in a screen. The relevant variables to be set are in the beginning of the bash script and they use the EVNT file from that is mentioned in the run_ana.py. 

relevant rivet analysis:
- rivet_anas/templates\*.cc
	- templatesnoptbins.cc may not work in the current state use templates_inclusive instead
- rivet_anas/angcoeff\*.cc

weird solutions to be aware of mostly because C++ doesnt support JSON:
- ptbins is saved as path/ptbins.txt where each line is a binedge in angcoeff\*.cc --> might need to check path
- ptbins file read by templates\*.cc --> need to check path
- templates/*.cc reads coeffs from a specified path !this will **not** be adjusted by the bash scripts --> need to check path!
	- coeffs.txt has the following form to be read correctly: each line = one coeff A_ij where i is the coeff number and j is the ptbinnumber and then the file has to be: 

A00
A01
A10
A11
A20
A21
A30
A31
A40
A41
A50
A51
A60
A61
A70
A71
A80
A81

A8 == xsec in the ptbin

when using the rivet_anas/read_angcoeffs.py this will be done correctly when using tprofiles called A_i for each coefficient.

- ptbins can be set in the angcoeff*.cc in the beginning
- the jsons have 3 keys --> bins (angular bins), yvals/zvals(1/2D histos), y/zerrs

u can always contact me with questions ... i know this isnt that much of a documentation but i hope this is a start 
